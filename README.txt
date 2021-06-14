How to do particle generation and detector simulation:

This tutorial is meant to be run on SLAC computing resources using the neutrino-jupyter/latest image. It can be run else where but an image with all the required packages would need to be configured. 

First navigate to the edep-sim_scripts folder. This is where the event generation and particle energy deposition simulation is done. It uses Geant4 macros and edep-sim. See the instructions in that folder. After running edep-sim, dumpTree.py must be run to convert to a .h5 file.

Second navigate to the larnd-sim_scripts folder. There are instructions there for how to run larnd-sim there. It takes as input the .h5 file produced at the last step. larnd-sim can be run in a jupyter notebook ('Pixel Induced Current.ipynb') or as a python script (simulate_pixels.py).
