This will describe how to generate particles using a Geant4 macro and edep-sim as a way to simulate the energy deposition in a detector.

'singlemuon.mac' is a simple Geant4 macro that produces a single monoenergetic muon from a point source at the center of the detector. This is just a simple example to use for testing. See 'macros' file for a few other examples.

We need to run edep-sim. In a terminal, run for example:

'edep-sim -u -e 10 -g Module0.gdml -o electron_test.root electron.mac'

'-e 10' generates 10 muons, change that as you'd like. '-g Module0.gdml' loads in the geometry for the detector we want to use.
'-o electron_test.root' specifies the name of the output file, should be a .root file. At the end the macro is specified.
'-u' does an update before running macros (I'm not sure how important this is, or what it does).

Once you produce the root file, it must be converted to a friendlier format. Run:

'python3 dumpTree.py electron_test.root electron_test.h5'

Which runs a python script that converts to an h5py format .h5. See 'explore_h5.ipynb' for an example of how to open and read h5py files. This is a crucial step (the converting) because larnd-sim, the rest of the detector simulation, uses this format as input. 
