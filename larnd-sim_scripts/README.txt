This describes how to run larnd-sim. This must be done after edep-sim is run, and a .h5 file has been created by dumpTree.py.

You need to git clone larnd-sim into this directory by running:

git clone https://github.com/DUNE/larnd-sim.git 

And you need to install it by running

python3 setup.py install

where you may need to include "--user" at the end if you do not have certain write permissions.

larnd-sim can be run two ways:

1) In the larnd-sim Examples folder, you can find Step-by-step simulation.ipynb. This is a way of running the software in a Jupyter notebook. This is also a great way to see what is going on more indepth and see some plots of the final products. 

2) Run the 'cli/simulate_pixels.py' script in a terminal. To see the usage, just type "python3 simulate_pixels.py".
Here is an example of running it:                                                    
                                                        
python3 cli/simulate_pixels.py \
--input_filename=input_file.h5 \
--detector_properties=larndsim/detector_properties/module0.yaml \
--pixel_layout=larndsim/pixel_layouts/multi_tile_layout-2.2.16.yaml \
--output_filename=output_file.h5 \
--response=response.npy

There are two yaml files required including one for pixel-level information and another that contains detector properties (both of which you may need to change for your detector).

See the slack channel @larnd-sim on DUNE for help. This example is for the Module0 prototype DUNE Near Detector. There also exists the corresponding files for a single TPC (HOWEVER, as of Nov 15, 2021, those files are not up to date. So someone will need to make them up to date...)

To run the event-display:

First convert the .h5 file produced in the previous step (not the one produced after the dumpTree step) to evd format:

'python3 event-display/to_evd_file.py --in_filename larndsim_test.h5 --out_filename larndsim_test_evd.h5 --geometry_file multi_tile_layout-2.2.16.yaml'

Then run the event-display:

'python3 event-display/quick_display.py -i larndsim_test_evd.h5'
