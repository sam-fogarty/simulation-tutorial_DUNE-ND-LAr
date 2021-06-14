import sys
import warnings
from datetime import datetime

import numpy as np
import fire
import h5py
import yaml

import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.colors import ListedColormap

from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

try:
    matplotlib.use("Qt5Agg")
except (ModuleNotFoundError, ImportError):
    matplotlib.use("TkAgg")
    print("Impossible to found Qt libraries, using Tkinter backend")

plt.ion()


class EventDisplay:

    def __init__(self, filename, geometry_file=None, nhits=1):
        f = h5py.File(filename, 'r')

        events = f['events']
        self.events = events[events['nhit'] > nhits]
        self.tracks = f['tracks'] if 'tracks' in f.keys() else None
        self.hits = f['hits']
        self.ext_trigs = f['ext_trigs'] if 'ext_trigs' in f.keys() else None
        self.info = f['info'].attrs

        self.fig = plt.figure(constrained_layout=False, figsize=(8.5, 6.))
        gs_xyzy = self.fig.add_gridspec(nrows=1, ncols=3, top=0.93, width_ratios=[1, 1, 0.05],
                                        left=0.15, right=0.5, bottom=0.58,
                                        hspace=0, wspace=0)

        ax_zy = self.fig.add_subplot(gs_xyzy[0])
        ax_xy = self.fig.add_subplot(gs_xyzy[1], sharey=ax_zy)

        cax = self.fig.add_subplot(gs_xyzy[2])
        ip = InsetPosition(ax_xy, [1.1, 0, 0.1, 1])
        cax.set_axes_locator(ip)
        self.cax = cax
        self.fig.text(
            0.04, 0.25, 'charge [10$^3$ e]', ha='center', va='center', rotation='vertical')

        gs_time = self.fig.add_gridspec(nrows=2, ncols=1,
                                        left=0.1, right=0.57, bottom=0.08, top=0.47,
                                        hspace=0.09)
        ax_time_1 = self.fig.add_subplot(gs_time[0])
        ax_time_2 = self.fig.add_subplot(gs_time[1], sharex=ax_time_1)

        gs_xyz = self.fig.add_gridspec(nrows=1, ncols=1,
                                       left=0.52, right=0.99, bottom=0.1, top=0.95,
                                       wspace=0)
        ax_xyz = self.fig.add_subplot(gs_xyz[0], projection='3d')
        ax_xyz.set_facecolor('none')

        if not geometry_file:
            geometry_file = self.info['geometry_file']

        with open(geometry_file, 'r') as gf:
            tile_layout = yaml.load(gf, Loader=yaml.FullLoader)

        mm2cm = 0.1
        pixel_pitch = tile_layout['pixel_pitch'] * mm2cm
        chip_channel_to_position = tile_layout['chip_channel_to_position']
        tile_chip_to_io = tile_layout['tile_chip_to_io']
        self.io_group_io_channel_to_tile = {}
        for tile in tile_chip_to_io:
            for chip in tile_chip_to_io[tile]:
                io_group_io_channel = tile_chip_to_io[tile][chip]
                io_group = io_group_io_channel//1000
                io_channel = io_group_io_channel % 1000
                self.io_group_io_channel_to_tile[(io_group, io_channel)] = tile

        cm2mm = 10

        xs = np.array(list(chip_channel_to_position.values()))[
            :, 0] * pixel_pitch * cm2mm
        ys = np.array(list(chip_channel_to_position.values()))[
            :, 1] * pixel_pitch * cm2mm
        tile_borders = np.zeros((2, 2))
        tpc_borders = np.zeros((0, 3, 2))
        tpc_centers = np.array(list(tile_layout['tpc_centers'].values()))
        tile_borders[0] = [-(max(xs)+pixel_pitch)/2, (max(xs)+pixel_pitch)/2]
        tile_borders[1] = [-(max(ys)+pixel_pitch)/2, (max(ys)+pixel_pitch)/2]

        tile_positions = np.array(list(tile_layout['tile_positions'].values()))
        tile_orientations = np.array(
            list(tile_layout['tile_orientations'].values()))
        self.tile_positions = tile_positions
        self.tile_orientations = tile_orientations
        tpcs = np.unique(tile_positions[:, 0])
        tpc_borders = np.zeros((len(tpcs), 3, 2))

        self.drift_length = abs(tile_positions[0][0])
        self.drift_time = self.drift_length / \
            self.info['vdrift']/self.info['clock_period']

        for itpc, tpc_id in enumerate(tpcs):
            this_tpc_tile = tile_positions[tile_positions[:, 0] == tpc_id]
            this_orientation = tile_orientations[tile_positions[:, 0] == tpc_id]

            x_border = min(this_tpc_tile[:, 2])+tile_borders[0][0]+tpc_centers[itpc][0], \
                max(this_tpc_tile[:, 2]) + \
                tile_borders[0][1]+tpc_centers[itpc][0]
            y_border = min(this_tpc_tile[:, 1])+tile_borders[1][0]+tpc_centers[itpc][1], \
                max(this_tpc_tile[:, 1]) + \
                tile_borders[1][1]+tpc_centers[itpc][1]
            z_border = min(this_tpc_tile[:, 0])+tpc_centers[itpc][2], \
                max(this_tpc_tile[:, 0])+self.drift_length * \
                this_orientation[:, 0][0]+tpc_centers[itpc][2]

            tpc_borders[itpc] = (x_border, y_border, z_border)

        self.tpc_borders = tpc_borders
        self.ax_xyz = ax_xyz
        self.ax_xy = ax_xy
        self.ax_time_1 = ax_time_1
        self.ax_time_2 = ax_time_2
        self.ax_zy = ax_zy

        self.run()

    def run(self):
        ev_id = 0

        while True:
            self.display_event(ev_id)

            user_input = input(
                'Next event (q to exit/enter for next/number to skip to position)?\n')
            if not user_input:
                ev_id += 1
            elif user_input[0].lower() == 'q':
                sys.exit()
            else:
                try:
                    ev_id = int(user_input)
                    print(ev_id)
                except ValueError:
                    print("Event number %s not valid" % user_input)
            if ev_id >= len(self.events):
                print("End of file")
                sys.exit()

    def _get_tile_id(self, io_group, io_channel):
        if (io_group, io_channel) in self.io_group_io_channel_to_tile:
            tile_id = self.io_group_io_channel_to_tile[io_group, io_channel]
        else:
            warnings.warn("IO group %i, IO channel %i not found" %
                          (io_group, io_channel))
            return 0

        return tile_id

    def _get_z_coordinate(self, io_group, io_channel, time):
        tile_id = self._get_tile_id(io_group, io_channel)

        z_anode = self.tile_positions[tile_id-1][0]
        drift_direction = self.tile_orientations[tile_id-1][0]

        return z_anode + time*self.info['vdrift']*self.info['clock_period']*drift_direction

    def get_event_start_time(self, event):
        """Estimate the event start time"""
        if event['n_ext_trigs']:
            # First Choice: Use earliest light system trigger in event
            return self.ext_trigs[event['ext_trig_ref']]['ts'][0]
        # Second Choice:
        #  Try to determine the start time from a 'bump' in charge.
        #  This is only valid if some part of the event hits one of the anodes.
        ticks_per_qsum = 10  # clock ticks per time bin
        t0_charge_threshold = 200.0  # Rough qsum threshold
        hit_ref = event['hit_ref']
        hits = self.hits[hit_ref]
        # determine charge vs time in enlarged window
        min_ts = np.amin(hits['ts'])
        max_ts = np.amax(hits['ts'])
        # If event long enough, calculate qsum vs time
        if (max_ts - min_ts) > ticks_per_qsum:
            time_bins = np.arange(min_ts-ticks_per_qsum,
                                  max_ts+ticks_per_qsum)
            # integrate q in sliding window to produce qsum profile
            #  histogram raw charge
            q_vs_t = np.histogram(hits['ts'],
                                  bins=time_bins,
                                  weights=hits['q'])[0]
            #  calculate rolling qsum
            qsum_vs_t = np.convolve(q_vs_t,
                                    np.ones(ticks_per_qsum, dtype=int),
                                    'valid')
            t0_bin_index = np.argmax(qsum_vs_t > t0_charge_threshold)
            t0_bin_index += ticks_per_qsum
            start_time = time_bins[t0_bin_index]
            # Check if qsum exceed threshold
            if start_time < max_ts:
                return start_time
        # Fallback is to use the first hit
        return event['ts_start']

    def set_axes(self):
        # self.ax_time_1.set_xticklabels([])
        # self.ax_time_1.set_xlim(0,self.drift_time)
        # self.ax_time_2.set_xlim(0,self.drift_time)
        self.ax_time_2.set_xlabel(r"timestamp [0.1 $\mathrm{\mu}$s]")
        self.ax_time_1.set_title("TPC 1", fontsize=10, x=0.5, y=0.75)
        self.ax_time_2.set_title("TPC 2", fontsize=10, x=0.5, y=0.75)

        self.ax_xy.set_xlim(np.min(self.tpc_borders[:, 2, :]), np.max(
            self.tpc_borders[:, 2, :]))
        self.ax_xy.set_ylim(np.min(self.tpc_borders[:, 1, :]), np.max(
            self.tpc_borders[:, 1, :]))
        self.ax_xy.set_aspect('equal')
        self.ax_xy.set_xlabel("x [mm]")
        for tk in self.ax_xy.get_yticklabels():
            tk.set_visible(False)
        # self.ax_xy.set_yticklabels([])

        self.ax_zy.set_xlim(np.min(self.tpc_borders[:, 0, :]), np.max(
            self.tpc_borders[:, 0, :]))
        self.ax_zy.set_ylim(np.min(self.tpc_borders[:, 1, :]), np.max(
            self.tpc_borders[:, 1, :]))
        self.ax_zy.set_aspect('equal')
        self.ax_zy.set_xlabel("z [mm]")
        self.ax_zy.set_ylabel("y [mm]")

        self.ax_zy.axvline(0, c='gray')

        anode1 = plt.Rectangle((self.tpc_borders[0][0][0], self.tpc_borders[0][1][0]),
                               self.tpc_borders[0][0][1] -
                               self.tpc_borders[0][0][0],
                               self.tpc_borders[0][1][1] -
                               self.tpc_borders[0][1][0],
                               linewidth=1, fc='none',
                               edgecolor='gray')
        self.ax_xyz.add_patch(anode1)
        art3d.pathpatch_2d_to_3d(anode1, z=self.tpc_borders[0][2][0], zdir="y")

        anode2 = plt.Rectangle((self.tpc_borders[0][0][0], self.tpc_borders[0][1][0]),
                               self.tpc_borders[0][0][1] -
                               self.tpc_borders[0][0][0],
                               self.tpc_borders[0][1][1] -
                               self.tpc_borders[0][1][0],
                               linewidth=1, fc='none',
                               edgecolor='gray')
        self.ax_xyz.add_patch(anode2)
        art3d.pathpatch_2d_to_3d(anode2, z=self.tpc_borders[1][2][0], zdir="y")

        cathode = plt.Rectangle((self.tpc_borders[0][0][0], self.tpc_borders[0][1][0]),
                                self.tpc_borders[0][0][1] -
                                self.tpc_borders[0][0][0],
                                self.tpc_borders[0][1][1] -
                                self.tpc_borders[0][1][0],
                                linewidth=1, fc='gray', alpha=0.25,
                                edgecolor='gray')
        self.ax_xyz.add_patch(cathode)
        art3d.pathpatch_2d_to_3d(cathode, z=0, zdir="y")

        self.ax_xyz.plot((self.tpc_borders[0][0][0], self.tpc_borders[0][0][0]),
                         (self.tpc_borders[0][2][0],
                          self.tpc_borders[1][2][0]),
                         (self.tpc_borders[0][1][0], self.tpc_borders[0][1][0]), lw=1, color='gray')

        self.ax_xyz.plot((self.tpc_borders[0][0][0], self.tpc_borders[0][0][0]),
                         (self.tpc_borders[0][2][0],
                          self.tpc_borders[1][2][0]),
                         (self.tpc_borders[0][1][1], self.tpc_borders[0][1][1]), lw=1, color='gray')

        self.ax_xyz.plot((self.tpc_borders[0][0][1], self.tpc_borders[0][0][1]),
                         (self.tpc_borders[0][2][0],
                          self.tpc_borders[1][2][0]),
                         (self.tpc_borders[0][1][0], self.tpc_borders[0][1][0]), lw=1, color='gray')

        self.ax_xyz.plot((self.tpc_borders[0][0][1], self.tpc_borders[0][0][1]),
                         (self.tpc_borders[0][2][0],
                          self.tpc_borders[1][2][0]),
                         (self.tpc_borders[0][1][1], self.tpc_borders[0][1][1]), lw=1, color='gray')

        self.ax_xyz.set_ylim(
            np.min(self.tpc_borders[:, 2, :]), np.max(self.tpc_borders[:, 2, :]))
        self.ax_xyz.set_xlim(
            np.min(self.tpc_borders[:, 0, :]), np.max(self.tpc_borders[:, 0, :]))
        self.ax_xyz.set_zlim(
            np.min(self.tpc_borders[:, 1, :]), np.max(self.tpc_borders[:, 1, :]))
        self.ax_xyz.grid(False)
        self.ax_xyz.set_xlabel("x [mm]")
        self.ax_xyz.set_ylabel("z [mm]")
        self.ax_xyz.set_zlabel("y [mm]")
        self.ax_xyz.set_box_aspect((2, 2, 4))
        self.ax_xyz.xaxis.set_major_locator(plt.MaxNLocator(3))
        self.ax_xyz.yaxis.set_major_locator(plt.MaxNLocator(3))
        self.ax_xyz.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        self.ax_xyz.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        self.ax_xyz.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        self.ax_xyz.zaxis.labelpad = 20

    def display_event(self, ev_id):
        self.ax_time_1.cla()
        self.ax_time_2.cla()
        self.ax_xyz.cla()
        self.ax_zy.cla()
        self.ax_xy.cla()
        self.cax.cla()
        self.set_axes()

        event = self.events[ev_id]
        event_datetime = datetime.utcfromtimestamp(
            event['unix_ts']).strftime('%Y-%m-%d %H:%M:%S')
        self.fig.suptitle("Event %i, ID %i - %s UTC" %
                          (ev_id, event['evid'], event_datetime))
        hit_ref = event['hit_ref']
        ext_trig_ref = event['ext_trig_ref']

        event_start_time = self.get_event_start_time(event)

        hits = self.hits[hit_ref]
        cmap = plt.cm.get_cmap('plasma')
        norm = matplotlib.colors.Normalize(
            vmin=min(self.hits[hit_ref]['q']), vmax=max(self.hits[hit_ref]['q']))
        mcharge = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        hits_anode1 = hits[hits['iogroup'] == 1]
        hits_anode2 = hits[hits['iogroup'] == 2]

        q_anode1 = hits_anode1['q'] * 0.250
        q_anode2 = hits_anode2['q'] * 0.250

        t_anode1 = hits_anode1['ts']-event_start_time
        t_anode2 = hits_anode2['ts']-event_start_time
        self.ax_time_1.hist(t_anode1, weights=q_anode1,
                            bins=200,  # np.linspace(0,self.drift_time,200),
                            histtype='step', label='binned')
        self.ax_time_2.hist(t_anode2, weights=q_anode2,
                            bins=200,  # np.linspace(0,self.drift_time,200),
                            histtype='step', label='binned')

        if q_anode1.any():
            self.ax_time_1.scatter(
                t_anode1, q_anode1, c='r', s=5, label='hits', zorder=999)
            self.ax_time_1.legend()

        if q_anode2.any():
            self.ax_time_2.scatter(
                t_anode2, q_anode2, c='r', s=5, label='hits', zorder=999)
            if not q_anode1.any():
                self.ax_time_2.legend()

        self.fig.colorbar(mcharge, cax=self.cax, label=r'Charge [$10^3$] e')

        if event['n_ext_trigs']:
            trig_delta = self.ext_trigs[ext_trig_ref]['ts']-event_start_time
            for trig in trig_delta:
                self.ax_time_1.axvline(x=trig, c='g')
                self.ax_time_2.axvline(x=trig, c='g')

        unassoc_hit_mask = np.ones(event['nhit']).astype(bool)

        if event['ntracks']:
            track_ref = event['track_ref']
            tracks = self.tracks[track_ref]
            track_start = tracks['start']
            track_end = tracks['end']
            for i, track in enumerate(tracks):

                hit_trk_ref = track['hit_ref']
                hits_trk = self.hits[hit_trk_ref]

                # Difference between the z coordinate using the event ts_start (used in the track fitter)
                # and the start time found by get_event_start_time
                z_correction = (self._get_z_coordinate(hits_trk['iogroup'][0], hits_trk['iochannel'][0], event_start_time)
                                - self._get_z_coordinate(hits_trk['iogroup'][0], hits_trk['iochannel'][0], event['ts_start']))

                self.ax_xy.plot((track_start[i][0], track_end[i][0]),
                                (track_start[i][1], track_end[i][1]),
                                c='C{}'.format(i+1), alpha=0.75, lw=1)

                self.ax_zy.plot((track_start[i][2], track_end[i][2]),
                                (track_start[i][1], track_end[i][1]),
                                c='C{}'.format(i+1), alpha=0.75, lw=1)

                hits_anode1 = hits_trk[hits_trk['iogroup'] == 1]
                hits_anode2 = hits_trk[hits_trk['iogroup'] == 2]

                self.ax_xy.scatter(hits_trk['px'], hits_trk['py'], lw=0.2, ec='C{}'.format(
                    i+1), c=cmap(norm(hits_trk['q'])), s=5, alpha=0.75)

                hitz = [self._get_z_coordinate(io_group, io_channel, time) for io_group, io_channel, time in zip(
                    hits_trk['iogroup'], hits_trk['iochannel'], hits_trk['ts']-track['t0'])]

                self.ax_zy.scatter(hitz, hits_trk['py'], lw=0.2, ec='C{}'.format(
                    i+1), c=cmap(norm(hits_trk['q'])), s=5, alpha=0.75)
                self.ax_xyz.scatter(hits_trk['px'], hitz, hits_trk['py'], lw=0.2, ec='C{}'.format(
                    i+1), c=cmap(norm(hits_trk['q'])), s=5, alpha=0.75)

                self.ax_xyz.plot((track_start[i][0], track_end[i][0]),
                                 (track_start[i][2], track_end[i][2]),
                                 (track_start[i][1], track_end[i][1]),
                                 c='C{}'.format(i+1), alpha=0.5, lw=4)

                unassoc_hit_mask[np.in1d(hits['hid'], hits_trk['hid'])] = 0

        if np.any(unassoc_hit_mask):
            unassoc_hits = hits[unassoc_hit_mask]
            BG = np.asarray([1., 1., 1., ])
            my_cmap = cmap(np.arange(cmap.N))
            alphas = np.linspace(0, 1, cmap.N)
            # Mix the colors with the background
            for i in range(cmap.N):
                my_cmap[i, :-1] = my_cmap[i, :-1] * \
                    alphas[i] + BG * (1.-alphas[i])
            my_cmap = ListedColormap(my_cmap)
            hitz = [self._get_z_coordinate(io_group, io_channel, time) for io_group, io_channel, time in zip(
                unassoc_hits['iogroup'], unassoc_hits['iochannel'], unassoc_hits['ts']-event_start_time)]
            self.ax_xyz.scatter(unassoc_hits['px'], hitz, unassoc_hits['py'], lw=0.2, ec='C0', c=my_cmap(
                norm(unassoc_hits['q'])), s=5, alpha=0.75)
            self.ax_xy.scatter(unassoc_hits['px'], unassoc_hits['py'], lw=0.2, ec='C0', c=my_cmap(
                norm(unassoc_hits['q'])), s=5, alpha=0.75)
            self.ax_zy.scatter(hitz, unassoc_hits['py'], lw=0.2, ec='C0', c=my_cmap(
                norm(unassoc_hits['q'])), s=5, alpha=0.75)


if __name__ == '__main__':
    fire.Fire(EventDisplay)
