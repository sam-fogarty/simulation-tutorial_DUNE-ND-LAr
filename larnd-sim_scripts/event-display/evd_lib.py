import h5py
import warnings
#import threading
import os
import numpy as np
import sklearn.cluster as cluster
import sklearn.decomposition as dcomp
from skimage.measure import LineModelND, ransac
from collections import defaultdict
import queue
import json
import time
import yaml

region_ref = h5py.special_dtype(ref=h5py.RegionReference)


class NonThread(object):
    def __init__(self, target):
        self.target = target

    def start(self):
        self.target()

    def join(self):
        pass

    def is_alive(self):
        return False

class ExternalTriggerFinder(object):
    '''
    A class to extract external triggers from packet arrays

    This class has two parameters: `pacman_trigger_enabled` and `larpix_trigger_channels`

    The parameter `pacman_trigger_enabled` configures the `ExternalTriggerFinder` to
    extract packets of `packet_type == 7` as external triggers

    The parameter `larpix_trigger_channels` configures the `ExternalTriggerFinder` to
    extract triggers on particular larpix channels as external triggers. To specify,
    this parameter should be a dict of `<chip-key>: [<channel id>]` pairs. A special
    chip key of `'All'` can be used in the event that all triggers on a particular
    channel of any chip key should be extracted as external triggers.

    You can access and set the parameters at initialization::

        etf = ExternalTriggerFinder(pacman_trigger_enabled=True, larpix_trigger_channels=dict())

    or via the getter/setters::

        etf.get_parameters() # dict(pacman_trigger_enabled=True, larpix_trigger_channels=dict())
        etf.get_parameters('pacman_trigger_enabled') # dict(pacman_trigger_enabled=True)

        etf.set_parameters(pacman_trigger_enabled=True, larpix_trigger_channels={'1-1-1':[0]})

    '''

    def __init__(self, pacman_trigger_enabled=True, larpix_trigger_channels=None):
        if larpix_trigger_channels is None:
            larpix_trigger_channels = dict()
        self._larpix_trigger_channels = larpix_trigger_channels
        self._pacman_trigger_enabled = pacman_trigger_enabled
        self.set_parameters()

    def get_parameters(self, *args):
        rv = dict()
        for key in ('pacman_trigger_enabled', 'larpix_trigger_channels'):
            if key in args or not args:
                rv[key] = getattr(self, '_{}'.format(key))
        return rv

    def set_parameters(self, **kwargs):
        self._pacman_trigger_enabled = kwargs.get(
            'pacman_trigger_enabled', self._pacman_trigger_enabled)
        self._larpix_trigger_channels = kwargs.get(
            'larpix_trigger_channels', self._larpix_trigger_channels)

    def fit(self, events, metadata=None):
        '''
        Pull external triggers from hit data within each event. No metadata is used.

        Trigger types are inherited from the pacman trigger type bits (with
        `pacman_trigger_enabled`) or are given a value of `-1` for larpix external triggers.

        :returns: a list of a list of dicts (one list for each event), each dict describes a single external trigger with the following keys: `ts`-trigger timestamp, `type`-trigger type, `mask`-mask for which packets within the event are included in the trigger

        '''
        if metadata is None:
            metadata = list()
        event_trigs = list()
        for event in events:
            event_trigs.append(list())

            if self._pacman_trigger_enabled:
                # first check for any pacman trigger packets
                trigger_mask = event['packet_type'] == 7
                trigger_mask[trigger_mask] = (
                    (event['trigger_type'][trigger_mask] >> 1) & 1).astype(bool)
                if np.any(trigger_mask) > 0:
                    for i, trigger in enumerate(event[trigger_mask]):
                        mask = np.zeros(len(event))
                        mask[i] = 1
                        event_trigs[-1].append(dict(
                            ts=trigger['timestamp'],
                            type=trigger['trigger_type'],
                            mask=mask.astype(bool)
                        ))

            if self._larpix_trigger_channels:
                # then check for larpix external triggers
                trigger_mask = np.zeros(len(event))
                for chip_key, channels in self._larpix_trigger_channels.items():
                    if chip_key == 'All':
                        key_mask = trigger_mask
                    else:
                        io_group, io_channel, chip_id = chip_key.split('-')
                        key_mask = np.logical_and.reduce((
                            io_group == event['io_group'],
                            io_channel == event['io_channel'],
                            chip_id == event['chip_id'],
                            trigger_mask
                        ))
                    for channel in channels:
                        trigger_mask = np.logical_or(np.logical_and(
                            event['channel_id'] == channel, key_mask), trigger_mask)
                trigger_mask = np.logical_and(
                    event['packet_type'] == 0, trigger_mask)
                if np.any(trigger_mask) > 0:
                    timestamps = event[trigger_mask]['timestamps']
                    split_indices = np.argwhere(
                        np.abs(np.diff(timestamps > 0)))
                    for idx, trigger in zip(split_indices, np.split(event[trigger_mask], split_indices)):
                        mask = np.zeros(len(event))
                        mask[trigger_mask][idx:idx+len(trigger)] = 1
                        event_trigs[-1].append(dict(
                            ts=np.mean(trigger['timestamp']),
                            type=-1,
                            mask=mask.astype(bool)
                        ))
        return event_trigs


class TrackFitter(object):
    '''
    A class to extract tracks from packet arrays

    You can access and set the parameters at initialization::

        tf = TrackFitter(vd=1.648, clock_period=0.1, ...)

    or via the getter/setters::

        tf.get_parameters() # dict(vd=1.648, clock_period=0.1, ...)
        tf.get_parameters('vd') # dict(vd=1.648)

        tf.set_parameters(vd=1.7)

    '''

    def __init__(self, dbscan_eps=14, dbscan_min_samples=5, vd=1.648, clock_period=0.1,
                 ransac_min_samples=2, ransac_residual_threshold=8, ransac_max_trials=100):
        self._vd = vd
        self._clock_period = clock_period
        self._dbscan_eps = dbscan_eps
        self._dbscan_min_samples = dbscan_min_samples
        self._ransac_min_samples = ransac_min_samples
        self._ransac_residual_threshold = ransac_residual_threshold
        self._ransac_max_trials = ransac_max_trials

        self.set_parameters()
        self.pca = dcomp.PCA(n_components=1)

    def get_parameters(self, *args):
        rv = dict()
        for key in ('vd', 'clock_period', 'z_scale',
                    'dbscan_eps', 'dbscan_min_samples',
                    'ransac_min_samples', 'ransac_residual_threshold', 'ransac_max_trials'):
            if key in args or not args:
                rv[key] = getattr(self, '_{}'.format(key))
        return rv

    def set_parameters(self, **kwargs):
        '''
        Sets parameters used in fitting tracks::

            vd:                        drift velocity [mm/us]
            clock_period:              clock period for timestamp [us]

            dbscan_eps:                epsilon used for clustering [mm]
            dbscan_min_samples:        min samples used for clustering

            ransac_min_samples:        min samples used for outlier detection
            ransac_residual_threshold: residual threshold used for outlier detection [mm]
            ransac_max_trials:         max trials used for outlier detection

        '''
        self._vd = kwargs.get('vd', self._vd)
        self._clock_period = kwargs.get('clock_period', self._clock_period)
        self._z_scale = self._vd * self._clock_period

        self._dbscan_eps = kwargs.get('dbscan_eps', self._dbscan_eps)
        self._dbscan_min_samples = kwargs.get(
            'dbscan_min_samples', self._dbscan_min_samples)
        self.dbscan = cluster.DBSCAN(
            eps=self._dbscan_eps, min_samples=self._dbscan_min_samples)

        self._min_samples = kwargs.get(
            'ransac_min_samples', self._ransac_min_samples)
        self._residual_threshold = kwargs.get(
            'ransac_residual_threshold', self._ransac_residual_threshold)
        self._max_trials = kwargs.get(
            'ransac_max_trials', self._ransac_max_trials)

    def _plot_dbscan(self, xyz):
        print('plotting for dbscan tuning...')
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        def refit(self, eps):
            self.dbscan.set_params(eps=eps)
            labels = self.dbscan.fit(xyz).labels_
            return len(np.unique(labels[labels != -1]))
        eps = np.linspace(0.1, 50, 100)
        n = [refit(self, e) for e in eps]
        plt.ion()
        fig = plt.figure()
        fig.add_subplot(1, 2, 1)
        plt.plot(eps, n)
        ax = fig.add_subplot(1, 2, 2, projection='3d')
        ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2])
        plt.show()
        plt.pause(5)

    def _do_dbscan(self, xyz, mask):
        clustering = self.dbscan.fit(xyz[mask])
        track_ids = np.zeros(len(mask))-1
        track_ids[mask] = clustering.labels_
        return track_ids

    def _do_ransac(self, xyz, mask):
        model_robust, inliers = ransac(xyz[mask], LineModelND,
                                       min_samples=self._ransac_min_samples,
                                       residual_threshold=self._ransac_residual_threshold,
                                       max_trials=self._ransac_max_trials)
        return inliers

    def _do_pca(self, xyz, mask):
        centroid = np.mean(xyz[mask], axis=0)
        pca = self.pca.fit(xyz[mask] - centroid)
        axis = pca.components_[0] / np.linalg.norm(pca.components_[0])
        return centroid, axis

    def _get_z_coordinate(self, io_to_tile, tile_geometry, io_group, io_channel, this_time):
        try:
            tile_id = io_to_tile[io_group, io_channel]
        except KeyError:
            return 0

        z_anode = tile_geometry[tile_id][0][0]
        drift_direction = tile_geometry[tile_id][1][0]

        return z_anode + this_time*self._z_scale*drift_direction

    def fit(self, event, metadata=None, plot=False):
        '''
        Extract tracks from a given event packet array

        Accepts geometry metadata and external trigger metadata (optional).
        Geometry should be specified as a dict of `(chip_id, channel_id)` pairs,
        and external triggers should be specified with a list of dicts containing
        `type` and `ts` keys.

        :returns: list of dicts (one for each track) containing keys: `track_id`-unique id within event, `mask`-mask for which packets are included in track, `centroid`-x,y,z centroid of track relative to `t0`, `axis`-track x,y,z axis, `residual`-x,y,z residuals, `length`-track length, `start`-x,y,z,t of track start point, `end`-x,y,z,t of track end point, `t0`-t0 timestamp used for track

        '''
        if metadata is None:
            metadata = dict()
        trigs = metadata.get('trigs', list())
        geometry = metadata.get('geometry', dict())

#         if trigs_list is None:
#             trigs_list = [list()]*len(events)
#         for event,trigs in zip(events,trigs_list):
#             event_tracks.append(list())
        if not len(geometry.keys()):
            return list()  # continue
        event_tracks = list()
        if len(event) < 2:
            return list()  # continue
        if trigs:
            t0 = np.min([trig['ts'] for trig in trigs]).astype(int)
            t0_type = 1
        else:
            t0 = event['timestamp'][0].astype(int)
            t0_type = 0
        iter_mask = np.ones(len(event)).astype(bool)
        while True:
            if metadata['tile_geometry']:
                xyz = np.array([(*geometry[(io_group, io_channel, chip_id, channel_id)],
                                 self._get_z_coordinate(metadata['io_to_tile'],
                                                        metadata['tile_geometry'],
                                                        io_group,
                                                        io_channel,
                                                        ts-t0))
                                for io_group, io_channel, chip_id, channel_id, ts in zip(event['io_group'],
                                                                                         event['io_channel'],
                                                                                         event['chip_id'],
                                                                                         event['channel_id'],
                                                                                         event['timestamp'].astype(int))])
            else:
                xyz = np.array([(*geometry[(1, 1, chip_id, channel_id)], (ts-t0)*self._z_scale) for chip_id,
                               channel_id, ts in zip(event['chip_id'], event['channel_id'], event['timestamp'].astype(int))])
            # dbscan to find clusters
            track_ids = self._do_dbscan(xyz, iter_mask)
            if plot:
                self._plot_dbscan(xyz)

            for track_id in np.unique(track_ids):
                if track_id == -1:
                    continue
                mask = track_ids == track_id
                if np.sum(mask) <= self._ransac_min_samples:
                    continue
                # ransac for linear hits
                inliers = self._do_ransac(xyz, mask)
                mask[mask] = inliers

                if np.sum(mask) < 2:
                    continue
                # PCA on central hits
                centroid, axis = self._do_pca(xyz, mask)
                r_min, r_max = self._projected_limits(
                    centroid, axis, xyz[mask])
                residual = self._track_residual(centroid, axis, xyz[mask])

                # convert back to global time
                r_min = np.append(r_min, r_min[-1]/self._z_scale)
                r_max = np.append(r_max, r_max[-1]/self._z_scale)

#                 event_tracks[-1].append(dict(
                event_tracks.append(dict(
                    track_id=track_id,
                    mask=mask,
                    centroid=centroid,
                    axis=axis,
                    residual=residual,
                    length=np.linalg.norm(r_max[:3]-r_min[:3]),
                    start=r_min,
                    end=r_max,
                    t0=t0,
                    t0_type=t0_type
                ))
                iter_mask[mask] = 0

            if np.all(track_ids == -1) or not np.any(iter_mask):
                break
        return event_tracks

    def _projected_limits(self, centroid, axis, xyz):
        s = np.dot((xyz - centroid), axis)
        xyz_min, xyz_max = np.amin(xyz, axis=0), np.amax(xyz, axis=0)
        r_max = np.clip(centroid + axis * np.max(s), xyz_min, xyz_max)
        r_min = np.clip(centroid + axis * np.min(s), xyz_min, xyz_max)
        return r_min, r_max

    def _track_residual(self, centroid, axis, xyz):
        s = np.dot((xyz - centroid), axis)
        res = np.abs(xyz - (centroid + np.outer(s, axis)))
        return np.mean(res, axis=0)

    def theta(self, axis):
        return np.arctan2(np.linalg.norm(axis[:2]), axis[-1])

    def phi(self, axis):
        return np.arctan2(axis[1], axis[0])

    def xyp(self, axis, centroid):
        if axis[-1] == 0:
            return centroid[:2]
        s = -centroid[-1] / axis[-1]
        return (centroid + axis * s)[:2]


class LArPixEVDFile(object):
    dtype_desc = {
        'info': None,
        'hits': [
            ('hid', 'i8'),
            ('px', 'f8'), ('py', 'f8'), ('ts', 'i8'), ('q', 'f8'),
            ('iochannel', 'i8'), ('iogroup',
                                  'i8'), ('chipid', 'i8'), ('channelid', 'i8'),
            ('geom', 'i8'), ('event_ref', region_ref), ('q_raw', 'f8')
        ],
        'events': [
            ('evid', 'i8'), ('track_ref', region_ref), ('hit_ref', region_ref),
            ('nhit', 'i8'), ('q', 'f8'), ('ts_start', 'i8'), ('ts_end', 'i8'),
            ('ntracks', 'i8'), ('ext_trig_ref', region_ref), ('n_ext_trigs', 'i8'),
            ('unix_ts', 'u8'), ('q_raw', 'f8')
        ],
        'tracks': [
            ('track_id', 'i8'), ('event_ref', region_ref), ('hit_ref', region_ref),
            ('theta', 'f8'), ('t0', 'i8'),
            ('phi', 'f8'), ('xp', 'f8'), ('yp', 'f8'), ('nhit', 'i8'),
            ('q', 'f8'), ('ts_start', 'i8'), ('ts_end', 'i8'),
            ('residual', 'f8', (3,)), ('length', 'f8'), ('start', 'f8', (4,)),
            ('end', 'f8', (4,)), ('q_raw', 'f8'), ('t0_type', 'u1')
        ],
        'ext_trigs': [
            ('trig_id', 'i8'), ('event_ref',
                                region_ref), ('ts', 'i8'), ('type', 'u1')
        ]
    }

    @staticmethod
    def _default_pxy():
        return (0., 0.)

    @staticmethod
    def _rotate_pixel(pixel_pos, tile_orientation):
        return pixel_pos[0]*tile_orientation[2], pixel_pos[1]*tile_orientation[1]

    def __init__(self, filename, source_file=None, configuration_file=None, geometry_file=None,
                 pedestal_file=None, builder_config=None, fitter_config=None, buffer_len=2048,
                 verbose=False, fit_tracks=True, trigger_finder_config=None, find_triggers=True,
                 cores=1, force=False, electron_lifetime_file=None):
        self.verbose = verbose
        self.is_open = True
        if os.path.exists(filename) and not force:
            raise OSError('{} exists!'.format(filename))
        elif os.path.exists(filename):
            print('deleting existing file {}...'.format(filename))
            os.remove(filename)
        self.h5_filename = filename
        self.out_buffer = list()
        self.buffer_len = buffer_len
        self.fit_tracks = fit_tracks
        self.cores = cores
        self.is_multi_tile = False

        # geometry lookup
        self.geometry = defaultdict(self._default_pxy)
        self.geometry_file = geometry_file
        self.tile_geometry = defaultdict(int)
        if geometry_file is not None:
            with open(geometry_file) as gf:
                geometry_yaml = yaml.load(gf, Loader=yaml.FullLoader)

            if 'multitile_layout_version' in geometry_yaml.keys():
                pixel_pitch = geometry_yaml['pixel_pitch']
                self.is_multi_tile = True
                chip_channel_to_position = geometry_yaml['chip_channel_to_position']
                tile_orientations = geometry_yaml['tile_orientations']
                tile_positions = geometry_yaml['tile_positions']
                tpc_centers = geometry_yaml['tpc_centers']
                tile_indeces = geometry_yaml['tile_indeces']
                xs = np.array(list(chip_channel_to_position.values()))[
                    :, 0] * pixel_pitch
                ys = np.array(list(chip_channel_to_position.values()))[
                    :, 1] * pixel_pitch
                x_size = max(xs)-min(xs)+pixel_pitch
                y_size = max(ys)-min(ys)+pixel_pitch

                self.io_group_io_channel_to_tile = {}
                for tile in geometry_yaml['tile_chip_to_io']:
                    tile_orientation = tile_orientations[tile]
                    self.tile_geometry[tile] = tile_positions[tile], tile_orientations[tile]
                    for chip in geometry_yaml['tile_chip_to_io'][tile]:
                        io_group_io_channel = geometry_yaml['tile_chip_to_io'][tile][chip]
                        io_group = io_group_io_channel//1000
                        io_channel = io_group_io_channel % 1000
                        self.io_group_io_channel_to_tile[(
                            io_group, io_channel)] = tile

                    for chip_channel in geometry_yaml['chip_channel_to_position']:
                        chip = chip_channel // 1000
                        channel = chip_channel % 1000
                        try:
                            io_group_io_channel = geometry_yaml['tile_chip_to_io'][tile][chip]
                        except KeyError:
                            #print("Chip %i on tile %i not present in network" % (chip,tile))
                            continue

                        io_group = io_group_io_channel // 1000
                        io_channel = io_group_io_channel % 1000
                        x = chip_channel_to_position[chip_channel][0] * \
                            pixel_pitch + pixel_pitch / 2 - x_size / 2
                        y = chip_channel_to_position[chip_channel][1] * \
                            pixel_pitch + pixel_pitch / 2 - y_size / 2

                        x, y = self._rotate_pixel((x, y), tile_orientation)
                        x += tile_positions[tile][2] + \
                            tpc_centers[tile_indeces[tile][1]][0]
                        y += tile_positions[tile][1] + \
                            tpc_centers[tile_indeces[tile][1]][1]

                        self.geometry[(io_group, io_channel,
                                       chip, channel)] = x, y

            else:
                import larpixgeometry.layouts
                geo = larpixgeometry.layouts.load(
                    self.geometry_file)  # open geometry yaml file
                for chip, pixels in geo['chips']:
                    for channel, pixel_id in enumerate(pixels):
                        if pixel_id is not None:
                            self.geometry[(1, 1, chip, channel)
                                          ] = geo['pixels'][pixel_id][1:3]

        # pedestal lookup
        self.pedestal = defaultdict(lambda: dict(
            pedestal_mv=580
        ))
        self.pedestal_file = pedestal_file
        if pedestal_file is not None:
            with open(self.pedestal_file, 'r') as infile:
                for key, value in json.load(infile).items():
                    self.pedestal[key] = value

        # chip configuration lookup
        self.configuration = defaultdict(lambda: dict(
            vref_mv=1300,
            vcm_mv=288
        ))
        self.configuration_file = configuration_file
        if configuration_file is not None:
            with open(self.configuration_file, 'r') as infile:
                for key, value in json.load(infile).items():
                    self.configuration[key] = value

        # electron lifetime lookup
        self.electron_lifetime_f = lambda unix, ts: 1.
        self.electron_lifetime_file = electron_lifetime_file
        if electron_lifetime_file is not None and fitter_config:
            if electron_lifetime_file[-5:] == '.root':
                import ROOT
                infile = ROOT.TFile(self.electron_lifetime_file, 'r')
                tg = infile.Get('lifetime')
                lt_unix = np.array([x for x in tg.GetX()])
                # convert from ms to clock ticks
                lt_tau = np.array([y for y in tg.GetY()]) * \
                    1e3 / (fitter_config.get('clock_period', 0.1))
                self.electron_lifetime_f = lambda unix, ts: np.exp(
                    ts.astype(float)/np.interp(unix.astype(float), lt_unix, lt_tau))
                infile.Close()
            else:
                self.electron_lifetime_file = None

        self.source_file = source_file

        fitter_config = fitter_config if fitter_config else dict()
        self.track_fitter = TrackFitter(**fitter_config)

        trigger_finder_config = trigger_finder_config if trigger_finder_config else dict()
        self.trigger_finder = ExternalTriggerFinder(**trigger_finder_config)

        self._queue = queue.Queue(2)
#         self._outfile_worker = multiprocessing.Process(target=self._parse_events_array)
#         self._outfile_worker = threading.Thread(target=self._parse_events_array)
        self._outfile_worker = NonThread(target=self._parse_events_array)

        self._create_datasets()
        self._write_metadata(dict(
            info=dict(
                source_file=self.source_file,
                vdrift=fitter_config['vd'] if fitter_config else 1.648,
                clock_period=fitter_config['clock_period'] if fitter_config else 0.1,
                configuration_file=self.configuration_file if self.configuration_file else '',
                pedestal_file=self.pedestal_file if self.pedestal_file else '',
                geometry_file=self.geometry_file if self.geometry_file else '',
                electron_lifetime_file=self.electron_lifetime_file if self.electron_lifetime_file else ''
            ),
            hits=dict(),
            events=builder_config if builder_config else dict(),
            tracks=self.track_fitter.get_parameters() if self.fit_tracks else dict()
        ))

    def append(self, event_array, unix_timestamp=0):
        self.out_buffer.append((event_array, unix_timestamp))
        if len(self.out_buffer) >= self.buffer_len:
            self.flush()

    def flush(self, block=True):
        if not len(self.out_buffer):
            return
        self._queue.put(self.out_buffer)
        if not self._outfile_worker.is_alive():
            #             self._outfile_worker = multiprocessing.Process(target=self._parse_events_array)
            #             self._outfile_worker = threading.Thread(target=self._parse_events_array)
            self._outfile_worker = NonThread(target=self._parse_events_array)
            self._outfile_worker.start()
        if block:
            self._outfile_worker.join()
        self.out_buffer = list()

    def close(self):
        self.flush(block=True)

    def _create_datasets(self):
        with h5py.File(self.h5_filename, 'a') as h5_file:
            for dataset_name, dataset_dtype in self.dtype_desc.items():
                if not dataset_dtype is None:
                    h5_file.create_dataset(
                        dataset_name, (0,), maxshape=(None,), dtype=dataset_dtype)
                else:
                    h5_file.create_group(dataset_name)

    def _write_metadata(self, metadata):
        with h5py.File(self.h5_filename, 'a') as h5_file:
            for name in metadata:
                for attr, value in metadata[name].items():
                    h5_file[name].attrs[attr] = value

    def _get_pixel_id(self, chip_id, channel_id):
        for chip_info in self.geometry['chips']:
            if chip_info[0] == chip_id:
                return chip_info[1][channel_id]
        return -1

    def _parse_events_array(self):
        last = time.time()
        with h5py.File(self.h5_filename, 'a') as h5_file:
            while not self._queue.empty():
                data_list = list()
                try:
                    data_list = self._queue.get(timeout=1)
                except queue.Empty:
                    break
                if self.verbose and time.time() > last+1:
                    print('{} chunks remaining...\r'.format(
                        self._queue.qsize()), end='')
                    last = time.time()

                events_list, unix_timestamps = zip(*data_list)

                # remove external triggers from events
                trigs_list = self.trigger_finder.fit(events_list)
                hit_masks = [~np.logical_or.reduce([trig['mask'] for trig in trigs])
                             if len(trigs) else np.ones(len(event)).astype(bool)
                             for event, trigs in zip(events_list, trigs_list)]
                events_list = [event[hit_mask]
                               for event, hit_mask in zip(events_list, hit_masks)]

                # run track fitting
#                 with multiprocessing.Pool(self.cores) as p:
    #                     tracks_list = self.track_fitter.fit(events_list, self.geometry, trigs_list=trigs_list) \
    #                         if self.fit_tracks else [list() for i in range(len(events_list))]
                tracks_list = [self.track_fitter.fit(
                    ev, dict(geometry=self.geometry, trigs=trigs_list[i], io_to_tile=self.io_group_io_channel_to_tile, tile_geometry=self.tile_geometry))
                    for i, ev in enumerate(events_list)] \
                    if self.fit_tracks else [list() for i in range(len(events_list))]

                # resize datasets
                events_dset = h5_file['events']
                tracks_dset = h5_file['tracks']
                hits_dset = h5_file['hits']
                ext_trigs_dset = h5_file['ext_trigs']
                events_idx = len(events_dset)
                tracks_idx = len(tracks_dset)
                hits_idx = len(hits_dset)
                ext_trigs_idx = len(ext_trigs_dset)
                events_dset.resize(
                    events_dset.shape[0] + len(events_list), axis=0)
                tracks_dset.resize(
                    tracks_dset.shape[0] + np.sum([len(tracks) for tracks in tracks_list]), axis=0)
                hits_dset.resize(
                    hits_dset.shape[0] + np.sum([len(event) for event in events_list]), axis=0)
                ext_trigs_dset.resize(
                    ext_trigs_dset.shape[0] + np.sum([len(trigs) for trigs in trigs_list]), axis=0)

                for unix_ts, event, tracks, trigs in zip(unix_timestamps, events_list, tracks_list, trigs_list):
                    events_dict = defaultdict(int)
                    tracks_dict = defaultdict(int)
                    hits_dict = defaultdict(int)
                    ext_trigs_dict = defaultdict(int)
                    # create event info
                    events_dict['unix_ts'] = np.array([unix_ts])
                    events_dict['nhit'] = np.array([len(event)])
                    events_dict['ntracks'] = np.array([len(tracks)])
                    events_dict['n_ext_trigs'] = np.array([len(trigs)])
                    events_dict['evid'] = np.array([events_idx])
                    events_dict['q'] = np.array([0.])
                    events_dict['q_raw'] = np.array([0.])
                    if len(event) and len(trigs):
                        events_dict['ts_start'] = np.array(
                            [min(event[0]['timestamp'], trigs[0]['ts'])])
                        events_dict['ts_end'] = np.array(
                            [max(event[-1]['timestamp'], trigs[-1]['ts'])])
                    elif len(event):
                        events_dict['ts_start'] = np.array(
                            [event[0]['timestamp']])
                        events_dict['ts_end'] = np.array(
                            [event[-1]['timestamp']])
                    elif len(trigs):
                        events_dict['ts_start'] = np.array([trigs[0]['ts']])
                        events_dict['ts_end'] = np.array([trigs[-1]['ts']])
                    # create trig info
                    if len(trigs):
                        ext_trigs_dict['trig_id'] = ext_trigs_idx + \
                            np.arange(len(trigs))
                        ext_trigs_dict['ts'] = np.array(
                            [trig['ts'] for trig in trigs])
                        ext_trigs_dict['type'] = np.array(
                            [trig['type'] for trig in trigs])
                    # create hit info
                    if len(event):
                        hits_dict['hid'] = hits_idx + np.arange(len(event))
                        hits_dict['ts'] = event['timestamp']
                        hits_dict['iogroup'] = event['io_group']
                        hits_dict['iochannel'] = event['io_channel']
                        hits_dict['chipid'] = event['chip_id']
                        hits_dict['channelid'] = event['channel_id']
                        hits_dict['geom'] = np.zeros(len(event))
                        hit_uniqueid = (((event['io_group'].astype(int))*256
                                         + event['io_channel'].astype(int))*256
                                        + event['chip_id'].astype(int))*64 \
                            + event['channel_id'].astype(int)
                        hit_uniqueid_str = hit_uniqueid.astype(str)
                        if self.is_multi_tile:
                            xy = np.array([self.geometry[(io_group, io_channel, chip_id, channel_id)]
                                           for io_group, io_channel, chip_id, channel_id in zip(event['io_group'], event['io_channel'], event['chip_id'], event['channel_id'])])
                        else:
                            xy = np.array([self.geometry[(
                                1, 1, (unique_id//64) % 256, unique_id % 64)] for unique_id in hit_uniqueid])

                        vref = np.array(
                            [self.configuration[unique_id]['vref_mv'] for unique_id in hit_uniqueid_str])
                        vcm = np.array([self.configuration[unique_id]['vcm_mv']
                                       for unique_id in hit_uniqueid_str])
                        ped = np.array([self.pedestal[unique_id]['pedestal_mv']
                                       for unique_id in hit_uniqueid_str])
                        hits_dict['px'] = xy[:, 0]
                        hits_dict['py'] = xy[:, 1]
                        q_raw = event['dataword']/256. * (vref-vcm) + vcm - ped
                        hits_dict['q_raw'] = q_raw
                        hits_dict['q'] = q_raw
                    # create track info
                    if len(tracks):
                        tracks_dict['track_id'] = tracks_idx + \
                            np.arange(len(tracks)).astype(int)
                        tracks_dict['nhit'] = np.zeros((len(tracks),))
                        tracks_dict['q_raw'] = np.zeros((len(tracks),))
                        tracks_dict['q'] = np.zeros((len(tracks),))
                        tracks_dict['ts_start'] = np.zeros((len(tracks),))
                        tracks_dict['ts_end'] = np.zeros((len(tracks),))
                        tracks_dict['theta'] = np.zeros((len(tracks),))
                        tracks_dict['phi'] = np.zeros((len(tracks),))
                        tracks_dict['xp'] = np.zeros((len(tracks),))
                        tracks_dict['yp'] = np.zeros((len(tracks),))
                        tracks_dict['residual'] = np.zeros((len(tracks), 3))
                        tracks_dict['length'] = np.zeros((len(tracks),))
                        tracks_dict['start'] = np.zeros((len(tracks), 4))
                        tracks_dict['end'] = np.zeros((len(tracks), 4))
                        tracks_dict['t0'] = np.zeros((len(tracks),))
                        tracks_dict['t0_type'] = np.zeros((len(tracks),))
                        for i, track in enumerate(tracks):
                            tracks_dict['t0'][i] = track['t0']
                            tracks_dict['t0_type'][i] = track['t0_type']
                            tracks_dict['nhit'][i] = np.sum(track['mask'])
                            tracks_dict['q_raw'][i] = np.sum(
                                hits_dict['q_raw'][track['mask']])
                            hits_dict['q'][track['mask']] = hits_dict['q_raw'][track['mask']] \
                                * self.electron_lifetime_f(
                                    events_dict['unix_ts'],
                                    hits_dict['ts'][track['mask']] -
                                track['t0']
                            )
                            tracks_dict['q'][i] = np.sum(
                                hits_dict['q'][track['mask']])
                            tracks_dict['theta'][i] = self.track_fitter.theta(
                                track['axis'])
                            tracks_dict['phi'][i] = self.track_fitter.phi(
                                track['axis'])
                            xp, yp = self.track_fitter.xyp(
                                track['axis'], track['centroid'])
                            tracks_dict['xp'][i] = xp
                            tracks_dict['yp'][i] = yp
                            tracks_dict['start'][i] = track['start']
                            tracks_dict['end'][i] = track['end']
                            tracks_dict['residual'][i] = track['residual']
                            tracks_dict['length'][i] = track['length']
                            tracks_dict['ts_start'][i] = hits_dict['ts'][track['mask']][0]
                            tracks_dict['ts_end'][i] = hits_dict['ts'][track['mask']][-1]

                    # calculate hit level info for events
                    if len(event):
                        events_dict['q_raw'] = np.array(
                            [np.sum(hits_dict['q_raw'])])
                        events_dict['q'] = np.array([np.sum(hits_dict['q'])])

                    # create references
                    event_ref = events_dset.regionref[events_idx]
                    track_ref = tracks_dset.regionref[tracks_idx:tracks_idx+len(
                        tracks)]
                    hit_ref = hits_dset.regionref[hits_idx:hits_idx+len(event)]
                    ext_trig_ref = ext_trigs_dset.regionref[ext_trigs_idx:ext_trigs_idx+len(
                        trigs)]

                    events_dict['track_ref'] = np.array([track_ref])
                    events_dict['hit_ref'] = np.array([hit_ref])
                    events_dict['ext_trig_ref'] = np.array([ext_trig_ref])
                    if len(trigs):
                        ext_trigs_dict['event_ref'] = np.array(
                            [event_ref]*len(trigs))
                    if len(tracks):
                        tracks_dict['event_ref'] = np.array(
                            [event_ref]*len(tracks))
                        tracks_dict['hit_ref'] = np.array(
                            [hits_dset.regionref[hits_dict['hid'][track['mask']]] for track in tracks])
                    if len(event):
                        hits_dict['event_ref'] = np.array(
                            [event_ref]*len(event))

                    self._fill(events_dset, events_idx, 1, **events_dict)
                    self._fill(ext_trigs_dset, ext_trigs_idx,
                               len(trigs), **ext_trigs_dict)
                    self._fill(tracks_dset, tracks_idx,
                               len(tracks), **tracks_dict)
                    self._fill(hits_dset, hits_idx, len(event), **hits_dict)

                    events_idx += 1
                    ext_trigs_idx += len(trigs)
                    tracks_idx += len(tracks)
                    hits_idx += len(event)

                self._queue.task_done()
        if self.verbose:
            print('0 chunks remaining...')

    def _fill(self, dset, start_idx, n, **kwargs):
        if n == 0:
            return
        data = np.zeros(n, dtype=self.dtype_desc[dset.name.strip('/')])
        for key, val in kwargs.items():
            data[key] = val
        dset[start_idx:start_idx+n] = data


def load_larpix_logfile(filename):
    datafile = h5py.File(filename, 'r')
    return datafile
