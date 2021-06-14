import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import numpy as np
import h5py
import argparse
import matplotlib
import sys
plt.ion()

def vol3d(x,y,z,q,*geom,name=None,fig=None,points=False):
    xyz = np.array(list(zip(x,y,z)))
    q = q+1e-9
    if not points:
        vox_q, edges = np.histogramdd(xyz, weights=q,
            bins=(
                np.linspace(geom[0],geom[1],
                    int((geom[1]-geom[0])/geom[-2])+1),
                np.linspace(geom[2],geom[3],
                    int((geom[3]-geom[2])/geom[-2])+1),
                np.linspace(geom[4],geom[5],
                    int((geom[5]-geom[4])/geom[-1])+1),
            ))
    norm = lambda x: np.clip((x - max(np.min(x),0.001)) / (np.max(x) - max(np.min(x),0.001)),0,1)
    cmap = plt.cm.get_cmap('plasma')
    if not points:
        vox_color = cmap(norm(vox_q))
        vox_color[..., 3] = norm(vox_q)
    else:
        vox_color = cmap(norm(q))
        vox_color[..., 3] = norm(q)

    ax = fig.add_subplot(1,2,2, projection='3d')
    if not points:
        ax.voxels(*np.meshgrid(*edges, indexing='ij'), vox_q, facecolors=vox_color)
    else:
        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],c=vox_color,alpha=0.5)
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('t [0.1us]')
    plt.xlim(geom[0],geom[1])
    plt.ylim(geom[2],geom[3])
    plt.tight_layout()

    plt.draw()
    return fig

def line3d(x,y,z,*geom,name=None,fig=None,points=False,c=None,edgecolors=None):
    xyz = np.array(list(zip(x,y,z)))

    if fig:
        axes = fig.get_axes()
        ax = axes[-1] if axes else fig.add_subplot(1,2,2, projection='3d')
    else:
        ax = fig.add_subplot(1,2,2, projection='3d')
    if not points:
        ax.plot(x,y,z,alpha=1,color=c)
    else:
        norm = lambda x: np.clip((x - max(np.min(x),0.001)) / (np.max(x) - max(np.min(x),0.001)),0,1)
        cmap = plt.cm.get_cmap('plasma')
        edgecolors = cmap(norm(edgecolors))
        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],s=5,alpha=0.75,c=c,edgecolors=edgecolors)
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('t [0.1us]')
    plt.xlim(geom[0],geom[1])
    plt.ylim(geom[2],geom[3])
    ax.set_zlim(geom[4],geom[5])
    plt.tight_layout()

    plt.draw()
    return fig

def plane3d(t,*geom,name=None,fig=None,alpha=0.2):
    if fig:
        axes = fig.get_axes()
        ax = axes[-1] if axes else fig.add_subplot(1,2,2, projection='3d')
    else:
        ax = fig.add_subplot(1,2,2, projection='3d')
    vtxs = [[
        (geom[0],geom[2]),
        (geom[0],geom[3]),
        (geom[1],geom[3]),
        (geom[1],geom[2])
    ]]

    p = PolyCollection(np.array(vtxs), alpha=alpha, facecolors='g')
    ax.add_collection3d(p,zs=t,zdir='z')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('t [0.1us]')

    plt.tight_layout()
    plt.draw()
    return fig

def line2d(x,y,*geom,name=None,fig=None,color=None):
    xy = np.array(list(zip(x,y)))

    ax = fig.gca()
    ax.plot(x,y,color=color)
    plt.tight_layout()

    plt.draw()
    return fig

def proj2d(x,y,q,*geom,name=None,fig=None):
    ax = fig.add_subplot(2,2,1)
    q = q+1e-9
    h = ax.hist2d(x,y,bins=(
        np.linspace(geom[0],geom[1],int((geom[1]-geom[0])/geom[-2])+1),
        np.linspace(geom[2],geom[3],int((geom[1]-geom[0])/geom[-2])+1)
        ),
        weights=q,
        cmin=0.0001,
        cmap='plasma'
    )
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.tight_layout()

    plt.draw()
    plt.colorbar(h[3],label='charge [ke]')
    return fig

def proj_time(t,q,*geom,name=None,fig=None):
    ax = fig.add_subplot(2,2,3)
    q = q+1e-9
    ax.hist(t, weights=q,
        bins=np.linspace(geom[4],geom[5],
                int((geom[5]-geom[4])/geom[-1])+1),
        histtype='step', label='binned')
    plt.xlabel('timestamp [0.1us]')
    plt.ylabel('charge [ke]')
    plt.tight_layout()

    plt.draw()
    return fig

def hit_times(t,q,*geom,name=None,fig=None):
    ax = fig.add_subplot(2,2,3)
    q = q+1e-9
    t,q = zip(*sorted(zip(t[t<geom[5]],q[t<geom[5]])))
    ax.plot(t,q,'r.', label='hits')
    plt.xlabel('timestamp [0.1us]')
    plt.ylabel('charge [ke]')
    plt.legend()
    plt.tight_layout()

    plt.draw()
    return fig

def trig_times(t,*geom,name=None,fig=None):
    ax = fig.add_subplot(2,2,3)
    ax.axvline(t,color='g',label='trigger')
    plt.xlabel('timestamp [0.1us]')
    plt.ylabel('charge [ke]')
    plt.legend()
    plt.tight_layout()

    plt.draw()
    return fig

def generate_plots(event, f, geom=[], fig=None, binned_3d=False):
    name = 'Event {}/{} ({})'.format(event['evid'],len(f['events']),f.filename)

    hits = f['hits']
    tracks = f['tracks'] if 'tracks' in f.keys() else None
    trigs = f['ext_trigs'] if 'ext_trigs' in f.keys() else None

    hit_ref = event['hit_ref']
    track_ref = event['track_ref'] if tracks else None
    ext_trig_ref = event['ext_trig_ref'] if trigs else None

    x = hits[hit_ref]['px']
    y = hits[hit_ref]['py']
    z = hits[hit_ref]['ts'] - event['ts_start']
    q = hits[hit_ref]['q'] * 0.250
    
    if not fig:
        fig = plt.figure(name,figsize=(8,6))

    if binned_3d:
        fig = vol3d(x,y,z,q,*geom,name=name,fig=fig)
    if tracks and 'ntracks' in event.dtype.names and event['ntracks']:
        track_start = tracks[track_ref]['start'][:,[0,1,3]]
        track_end = tracks[track_ref]['end'][:,[0,1,3]]
        track_start[:,2] += tracks[track_ref]['t0'] - event['ts_start']
        track_end[:,2] += tracks[track_ref]['t0'] - event['ts_start']
        unassoc_hit_mask = np.ones(event['nhit']).astype(bool)
        for i,track in enumerate(tracks[track_ref]):
            hit_ref = track['hit_ref']
            print('track {}\tlength: {:0.01f}mm\tevent frac: {:0.02f}\tnhit {}'.format(i,track['length'],track['nhit']/event['nhit'],track['nhit']))
            fig = line3d(hits[hit_ref]['px'], hits[hit_ref]['py'], hits[hit_ref]['ts']-event['ts_start'], *geom, edgecolors=hits[hit_ref]['q'], c='C{}'.format(i+1), name=name, fig=fig, points=True)
            fig = line3d((track_start[i,0],track_end[i,0]), (track_start[i,1],track_end[i,1]), (track_start[i,2],track_end[i,2]), *geom, name=name, fig=fig, c='C{}'.format(i+1))
            unassoc_hit_mask[np.in1d(hits[event['hit_ref']]['hid'],hits[hit_ref]['hid'])] = 0
        if np.sum(unassoc_hit_mask):
            unassoc_hits = hits[event['hit_ref']][unassoc_hit_mask]
            fig = line3d(unassoc_hits['px'], unassoc_hits['py'], unassoc_hits['ts']-event['ts_start'], *geom, edgecolors=unassoc_hits['q'], c='C0', name=name, fig=fig, points=True)
    else:
        unassoc_hits = hits[hit_ref]
        fig = line3d(unassoc_hits['px'], unassoc_hits['py'], unassoc_hits['ts']-event['ts_start'], *geom, edgecolors=unassoc_hits['q'], name=name, fig=fig, points=True, c='C0')
    if trigs and 'n_ext_trigs' in event.dtype.names and event['n_ext_trigs']:
        for ts in trigs[ext_trig_ref]['ts']:
            fig = plane3d(ts-event['ts_start'],*geom,name=name,fig=fig,alpha=0.2)

    fig = proj2d(x,y,q,*geom,name=name,fig=fig)
    if tracks and 'ntracks' in event.dtype.names and event['ntracks']:
        for i,se in enumerate(zip(track_start,track_end)):
            s,e = se
            fig = line2d((s[0],e[0]),(s[1],e[1]),*geom,name=name,fig=fig, color='C{}'.format(i+1))

    fig = proj_time(z,q,*geom,name=name,fig=fig)
    fig = hit_times(z,q,*geom,name=name,fig=fig)
    if trigs and 'n_ext_trigs' in event.dtype.names and event['n_ext_trigs']:
        for ts in trigs[ext_trig_ref]['ts']:
            fig = trig_times(ts-event['ts_start'],*geom,name=name,fig=fig)
    fig.canvas.set_window_title(name)
    return fig

def open_file(filename):
    return h5py.File(filename,'r')

def main(args):
    f = open_file(args.input)
    events = f['events']
    tracks = f['tracks'] if 'tracks' in f.keys() else None
    hits = f['hits']
    ext_trigs = f['ext_trigs'] if 'ext_trigs' in f.keys() else None
    fig = None
    ev = 0
    while True:
        print('displaying event {} with nhit_sel={}'.format(ev,args.nhit_sel))
        if ev >= np.sum(events['nhit'] > args.nhit_sel):
            sys.exit()
        event = events[events['nhit'] > args.nhit_sel][ev]
        print('Hits:',hits[event['hit_ref']])
        if tracks and 'ntracks' in event.dtype.names and event['ntracks']: print('Track:',tracks[event['track_ref']])
        if ext_trigs and 'n_ext_trigs' in event.dtype.names and event['n_ext_trigs']: print('Ext Triggger:',ext_trigs[event['ext_trig_ref']])
        print('Event:',event)
        fig = generate_plots(event, f, args.geom_limits, fig=fig, binned_3d=args.binned_3d)
        user_input = input('Next event (q to exit/enter for next/number to skip to position)?\n')
        print(user_input)
        if not len(user_input) or user_input[0] == '':
            ev += 1
        elif user_input[0].lower() == 'q':
            sys.exit()
        else:
            ev = int(user_input)
        plt.clf()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',required=True,help='''
    Input event display file
    ''')
    parser.add_argument('--nhit_sel',default=0, type=int, help='''
    Optional, sub-select on nhit greater than this value
    ''')
    parser.add_argument('--geom_limits', default=[-159.624,159.624,-159.624,159.624,0,1900,4.434,23], nargs=8, type=float, metavar=('XMIN','XMAX','YMIN','YMAX','TMIN','TMAX','PIXEL_PITCH','TIME_VOXEL'), help='''
    Optional, limits for geometry
    ''')
    parser.add_argument('--binned_3d',action='store_true',default=False,help='''
    Optional, use a 3D voxelization of charge in the 3D view (slower)
    ''')
    args = parser.parse_args()
    main(args)
