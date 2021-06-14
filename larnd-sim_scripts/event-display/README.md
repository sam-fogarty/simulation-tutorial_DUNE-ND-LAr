# Event display parser
This contains a simple library (`evd_lib.py`) and a script (`to_evd_file.py`) for converting v2.1 larpix-control hdf5 log files into a format that the larpix event display can interpret.

# Usage
```
python to_evd_file.py \
    --in_filename <input hdf5 log file, req.> \
    --out_filename <output hdf5 event file, req.> \
    --geometry_file <larpix geo yaml file for pixel location, opt.> \
    --pedestal_file <pedestal json file, opt.> \
    --configuration_file <configuration json file, opt.> \
    --buffer_size <packet buffer for sorting by timestamp, opt.> \
    --event_dt <max delta time for hits within an event [clk cycles], opt.> \
    --nhit_cut <min number of hits in an event, opt.>
```

## Input file
The input file must be a version 2.1 larpix-control hdf5 file (see `larpix-control/larpix/format/hdf5format` for a complete description)

## Output file
The output file is an hdf5 file with three datasets (`hits`,`events`,`tracks`). The `hits` dataset contains hit level information: pixel location, pixel id, charge, timestamp. The `tracks` dataset is largely unimplemented at this point, but represents single straight line segments within each event. The `events` dataset contains event level information where hits and tracks have been associated.

## Geometry file
This is a `larpix-geometry` generated yaml file that allows for association between the chip id / channel id to a unique pixel id and x,y location.

## Pedestal file
This is a json file formatted as::

    {
        <unique channel id>: {
            "pedestal_mv": <chip pedestal in mV>
        },
        ...
    }

The unique channel id is unique identifier for each channel, defined as `io_group * 4194304 + io_channel * 16384 + chip_id * 64 + channel_id`.

## Configuration file
This is a json file formatted as::

    {
        <unique chip id>: {
            "vref_mv": <vref of chip in mV>,
            "vcm_mv": <vcm of chip in mV>
        },
        ...
    }

The unique chip id is a unique identifier for each chip, defined as `io_group * 65536 + io_channel * 256 + chip_id`.

## Buffer size
To divide a continuous packet stream into discrete events required by the event display, the packet stream is by timestamp in small chunks and split whenever there is a large timestamp gap. This parameter sets the length of the chunk used to separate events. To avoid issues associated with the 512 bug of v2, this should generally be set >1536.

## Event delta time
A new event is created whenever there is a delta in two sequential timestamps greater than `event_dt` clk cycles.

## NHit cut
To speed up processing and reduce the output file size, a cut on the number of hits in an event can be placed. Any potential events that have few than `nhit_cut` triggers in them are excluded from the output data file.

