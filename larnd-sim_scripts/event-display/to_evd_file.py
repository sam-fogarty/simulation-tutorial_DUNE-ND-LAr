import numpy as np
import json
import argparse
import time

from evd_lib import *

_default_geometry_file         = None
_default_configuration_file    = None
_default_pedestal_file         = None
_default_buffer_size           = 38400
_default_event_dt              = 1820
_default_max_event_dt          = 3*1820
_default_nhit_cut              = 10
_default_max_packets           = -1
_default_dbscan_eps            = 14.
_default_dbscan_min_samples    = 5
_default_vd                    = 1.648
_default_clock_period          = 0.1
_default_skip_track_fit        = False
_default_external_trigger_conf = dict(
    pacman_trigger_enabled  = True,
    larpix_trigger_channels = {'All': [6]}
    )
_default_skip_trigger_finding   = False
_default_force                  = False
_default_electron_lifetime_file = None

def split_at_timestamp(timestamp,event,*args):
    '''
    Breaks event into two arrays at index where event['timestamp'] > timestamp
    Additional arrays can be specified with kwargs and will be split at the same
    index

    :returns: tuple of two event halves followed by any additional arrays (in pairs)
    '''
    args = list(args)
    timestamps = event['timestamp'].astype(int)
    indices = np.argwhere(timestamps > timestamp)
    if len(indices):
        idx = np.min(indices)
        args.insert(0,event)
        rv = [(arg[:idx], arg[idx:]) for arg in args]
        return tuple(v for vs in rv for v in vs)
    args.insert(0,event)
    rv = [(arg, np.array([], dtype=arg.dtype)) for arg in args]
    return tuple(v for vs in rv for v in vs)

def add_event(evd_file,event,unix_ts,max_event_dt=0):
    event,event_buffer,unix_ts,unix_ts_buffer = split_at_timestamp(
        min(event['timestamp'].astype(int)) + max_event_dt,
        event,
        unix_ts
    )
    event_timestamp,counts = np.unique(unix_ts['timestamp'], return_counts=True)
    evd_file.append(
        event,
        unix_timestamp=event_timestamp[np.argmax(counts)]
    )
    return event,event_buffer,unix_ts,unix_ts_buffer

def main(in_filename, out_filename, *args,
         configuration_file=_default_configuration_file, geometry_file=_default_geometry_file,
         pedestal_file=_default_pedestal_file,
         buffer_size=_default_buffer_size, event_dt=_default_event_dt, max_event_dt=_default_max_event_dt,
         nhit_cut=_default_nhit_cut,
         max_packets=_default_max_packets,
         dbscan_eps=_default_dbscan_eps, dbscan_min_samples=_default_dbscan_min_samples,
         vd=_default_vd, clock_period=_default_clock_period,
         skip_track_fit=_default_skip_track_fit,
         external_trigger_conf=_default_external_trigger_conf,
         skip_trigger_finding=_default_skip_trigger_finding,
         force=_default_force,
         electron_lifetime_file=_default_electron_lifetime_file,
         **kwargs):
    # load larpix file
    larpix_logfile = load_larpix_logfile(in_filename)
    packets        = larpix_logfile['packets']

    # create buffered output file
    evd_file = LArPixEVDFile(
        out_filename,
        source_file        = in_filename,
        geometry_file      = geometry_file,
        configuration_file = configuration_file,
        pedestal_file      = pedestal_file,
        builder_config     = dict(
            buffer_size = buffer_size,
            event_dt    = event_dt,
            nhit_cut    = nhit_cut,
            max_packets = max_packets
            ),
        fitter_config = dict(
            vd                 = vd,
            clock_period       = clock_period,
            dbscan_eps         = dbscan_eps,
            dbscan_min_samples = dbscan_min_samples
            ),
        fit_tracks             = not skip_track_fit,
        trigger_finder_config  = external_trigger_conf,
        find_triggers          = not skip_trigger_finding,
        force                  = force,
        electron_lifetime_file = electron_lifetime_file
        )
    event_counter  = 0
    packet_counter = 0

    # remove configuration messages and packets with bad parity
    good_parity_mask      = packets['valid_parity'][:]
    data_packet_mask      = packets['packet_type'][:] == 0
    trigger_packet_mask   = packets['packet_type'][:] == 7
    sync_packet_mask      = packets['packet_type'][:] == 6
    timestamp_packet_mask = packets['packet_type'][:] == 4
    mask = np.logical_and(good_parity_mask, data_packet_mask)
    mask = np.logical_or(mask,timestamp_packet_mask)
    if 'pacman_trigger_enabled' in external_trigger_conf and external_trigger_conf['pacman_trigger_enabled']:
        mask = np.logical_or(mask,trigger_packet_mask)
    del good_parity_mask
    del data_packet_mask
    del trigger_packet_mask
    del sync_packet_mask

    n_packets      = len(packets) #int(np.sum(mask))
    start_idx      = 0
    end_idx        = buffer_size
    event_buffer   = np.array([],dtype=packets.dtype)
    unix_ts_buffer = np.array([],dtype=packets.dtype)
    start_time     = time.time()
    last_unix_ts   = np.array(packets[timestamp_packet_mask][0],dtype=packets.dtype)
    del timestamp_packet_mask
    while start_idx < n_packets and (max_packets < 0 or start_idx < max_packets):
        # load a buffer of data
        block = packets[start_idx:min(end_idx,n_packets)]

        mask = (block['valid_parity'].astype(bool) & (block['packet_type'] == 0)) # data packets
        mask = mask | (block['packet_type'] == 4) # timestamp packets
        mask = mask | ((block['packet_type'] == 7) \
            & ('pacman_trigger_enabled' in external_trigger_conf) \
            & external_trigger_conf['pacman_trigger_enabled']) # external trigger packets

        packet_buffer = np.copy(block[mask])
        # packet_buffer = np.copy(packets[mask][start_idx:min(end_idx,n_packets)])
        packet_buffer = np.insert(packet_buffer, [0], last_unix_ts)

        # find unix timestamp groups
        ts_mask = packet_buffer['packet_type'] == 4
        ts_grps = np.split(packet_buffer, np.argwhere(ts_mask).flatten())
        unix_ts = np.concatenate([[ts_grp[0]]*len(ts_grp[1:]) for ts_grp in ts_grps if len(ts_grp) > 1], axis=0)
        packet_buffer = packet_buffer[~ts_mask]
        packet_buffer['timestamp'] = packet_buffer['timestamp'].astype(int) % (2**31) # ignore 32nd bit from pacman triggers
        last_unix_ts = unix_ts[-1]

        # sort packets to fix 512 bug
        sorted_idcs = np.argsort(packet_buffer,order='timestamp')
        packet_buffer = packet_buffer[sorted_idcs]
        unix_ts       = unix_ts[sorted_idcs]

        # cluster into events by delta t
        packet_dt = packet_buffer['timestamp'][1:] - packet_buffer['timestamp'][:-1]
        if len(event_buffer):
            packet_dt = np.insert(packet_dt, [0], packet_buffer['timestamp'][0] - event_buffer[-1]['timestamp'].astype(int))
        else:
            packet_dt = np.insert(packet_dt, [0], 0)
        event_idx     = np.argwhere(np.abs(packet_dt) > event_dt).flatten()
        events        = np.split(packet_buffer, event_idx)
        event_unix_ts = np.split(unix_ts, event_idx)

        for idx, event, unix_ts in zip(event_idx, events[:-1], event_unix_ts[:-1]):
            if idx == 0:
                while len(event_buffer) >= nhit_cut:
                    # current event buffer is a complete event
                    _,event_buffer,__,unix_ts_buffer = add_event(
                        evd_file,
                        event_buffer,
                        unix_ts_buffer,
                        max_event_dt=max_event_dt
                    )
                    event_counter  += 1
                event_buffer   = np.array([],dtype=packets.dtype)
                unix_ts_buffer = np.array([],dtype=packets.dtype)
            while len(event)+len(event_buffer) >= nhit_cut:
                if len(event_buffer):
                    # event found within packet buffer (combine with existing event buffer)
                    event   = np.append(event_buffer, event)
                    unix_ts = np.append(unix_ts_buffer, unix_ts)
                    event,event_buffer,unix_ts,unix_ts_buffer = add_event(
                        evd_file,
                        event,
                        unix_ts,
                        max_event_dt=max_event_dt
                    )
                    event_counter  += 1
                else:
                    event,event_buffer,unix_ts,unix_ts_buffer = add_event(
                        evd_file,
                        event,
                        unix_ts,
                        max_event_dt=max_event_dt
                    )
                    event_counter  += 1
                event   = np.array([],dtype=packets.dtype)
                unix_ts = np.array([],dtype=packets.dtype)
            event_buffer   = np.array([],dtype=packets.dtype)
            unix_ts_buffer = np.array([],dtype=packets.dtype)

        # keep any lingering packets for next iteration
        event_buffer   = np.array(events[-1],dtype=packets.dtype)
        unix_ts_buffer = np.array(event_unix_ts[-1],dtype=packets.dtype)

        # increment buffer indices
        start_idx = end_idx
        end_idx   = end_idx + buffer_size

        packet_counter += len(packet_buffer)
        print('packets parsed: {}/{}\tevents found: {}\ttime remaining {:0.02f}min...'.format(packet_counter, n_packets, event_counter, (time.time()-start_time)/packet_counter/60 * (n_packets-packet_counter)),end='\r')

    while len(event_buffer) >= nhit_cut:
        event,event_buffer,unix_ts,unix_ts_buffer = add_event(
            evd_file,
            event_buffer,
            unix_ts_buffer,
            max_event_dt=max_event_dt
        )
        event_counter  += 1
    print('packets parsed: {}/{}\tevents found: {}\ttime remaining {:0.02f}min...Done!'.format(packet_counter, n_packets, event_counter, (time.time()-start_time)/(packet_counter+1e-9)/60 * (n_packets-packet_counter)))

    print('finishing up...')
    evd_file.verbose = True
    evd_file.close()
    larpix_logfile.close()
    print('done')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_filename','-i',required=True,type=str)
    parser.add_argument('--out_filename','-o',required=True,type=str)
    parser.add_argument('--geometry_file','-g',default=_default_geometry_file,type=str,help='''default=%(default)s''')
    parser.add_argument('--pedestal_file','-p',default=_default_pedestal_file,type=str,help='''default=%(default)s''')
    parser.add_argument('--configuration_file','-c',default=_default_configuration_file,type=str,help='''default=%(default)s''')
    parser.add_argument('--buffer_size','-b',default=_default_buffer_size,type=int,help='''default=%(default)s''')
    parser.add_argument('--event_dt',default=_default_event_dt,type=int,help='''default=%(default)s''')
    parser.add_argument('--event_max_dt',default=_default_max_event_dt,type=int,help='''default=%(default)s''')
    parser.add_argument('--nhit_cut',default=_default_nhit_cut,type=int,help='''default=%(default)s''')
    parser.add_argument('--max_packets','-n',default=_default_max_packets,type=int,help='''default=%(default)s''')
    parser.add_argument('--vd',default=_default_vd,type=float,help='''default=%(default)s''')
    parser.add_argument('--clock_period',default=_default_clock_period,type=float,help='''default=%(default)s''')
    parser.add_argument('--dbscan_eps',default=_default_dbscan_eps,type=float,help='''default=%(default)s''')
    parser.add_argument('--dbscan_min_samples',default=_default_dbscan_min_samples,type=int,help='''default=%(default)s''')
    parser.add_argument('--skip_track_fit',action='store_true',help='''flag to skip track fitting''')
    parser.add_argument('--external_trigger_conf',default=_default_external_trigger_conf,type=json.loads,help='''config for external trigger finder, json-formatted dict with keys 'pacman_trigger_enabled' and 'larpix_trigger_channels', default=%(default)s''')
    parser.add_argument('--skip_trigger_finding',action='store_true',help='''flag to skip external trigger finding''')
    parser.add_argument('--force','-f',action='store_true',help='''overwrite file if it exists''')
    parser.add_argument('--electron_lifetime_file',default=_default_electron_lifetime_file,type=str,help='''file containing electron lifetime calibration''')
    args = parser.parse_args()
    main(**vars(args))
