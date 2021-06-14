This describes how to run larnd-sim. This must be done after edep-sim is run, and a .h5 file has been created by dumpTree.py. 

larnd-sim can be run two ways:

1) Run in a jupyter notebook; see 'Pixel Induced Current.ipynb'. This is great to see what is going on more indepth and see some plots of the final products. 

2) Run the 'simulate_pixels.py' script in a terminal. The usage is the following:

Usage: simulate_pixels.py INPUT_FILENAME PIXEL_LAYOUT DETECTOR_PROPERTIES <flags>
  optional flags:        --output_filename | --n_tracks
  
So, for example:
                            input edep-sim file                                       detector properties
                            |                                                         |
'python3 simulate_pixels.py ../edep-sim_scripts/test.h5 multi_tile_layout-2.2.16.yaml module0.yaml --output_filename larndsim_test.h5'
                                                        |                             
                                                        pixel layout for detector     
                             
The source for larnd-sim is https://github.com/DUNE/larnd-sim. See the slack channel @larnd-sim on DUNE for help. See the larndsim folder for the pixel layouts and detector properties. This example is for the Module0 prototype DUNE Near Detector. There also exists the corresponding files for a single TPC.

To run the event-display:

First convert the .h5 file produced in the previous step (not the one produced after the dumpTree step) to evd format:

'python3 event-display/to_evd_file.py --in_filename larndsim_test.h5 --out_filename larndsim_test_evd.h5 --geometry_file multi_tile_layout-2.2.16.yaml'

Then run the event-display:

'python3 event-display/quick_display.py -i larndsim_test_evd.h5'
