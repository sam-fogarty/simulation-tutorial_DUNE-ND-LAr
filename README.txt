How to do particle generation and detector simulation:

This tutorial can easily be run on SLAC computing resources using e.g. the neutrino-jupyter/latest image.  It can be run else where but an image with all the required packages would need to be configured. You will need to install edep-sim, ROOT, and all the packages larnd-sim requires. I recommend installing everything before attempting to run the tutorial.

First you should navigate to larnd-sim_scripts and pull and install larnd-sim from the git repository (info in a README there). This is to ensure you are running the latest version of larnd-sim, since it is still changing as it is adapted for the full ND-LAr. This software will work for any variation on ND-LAr, such as SingleCube and single Modules (e.g. Module 0). Support for the 2x2 and full ND-LAr is in-progress as of early November 2021 (see git repo for updates). As long as the simulation supports the detector of interest, you will just need the geometry for the detector, the pixel yaml file, and the detector properties yaml file. 

Now you can start simulating. Next, navigate to the edep-sim_scripts folder. This is where the event generation and particle energy deposition simulation is done. It uses Geant4 macros and edep-sim (but there are other options for event generation e.g. hepevt input, ROOTRACKER, etc). See the instructions in that folder, and I do highly recommend looking through the edep-sim git repository README for details on edep-sim. After running edep-sim, dumpTree.py must be run to convert the ROOT file to a .h5 file before larnd-sim can be run.

After you've made some edep-sim output files, navigate to the larnd-sim_scripts folder. There are instructions there for how to run larnd-sim there. It takes as input the .h5 file produced at the last step. larnd-sim can be run in a jupyter notebook or as a python script.

A side note: If you would like to improve or add on to larnd-sim (and are interested in applying those improvement/additions to the larnd-sim git repo), you need to keep in mind that only certain packages/commands will work with larnd-sim. This is because larnd-sim uses GPUs and the CUDA python package in order to speed up the software by a lot compared to just CPUs. As a result, only certain packages and python commands are compatible with CUDA. 

If you'd like to match charge and light events, either in simulation or data, who may want to look into module0_flow

https://github.com/peter-madigan/module0_flow
https://module0-flow.readthedocs.io/en/latest/

It is currently only working for module 0 (or possibly just a single module of ND-LAr). It "builds" charge events which basically means that it clusters and collects together charge packets into events. This is also a step that you could do on your own without this software, but module0_flow provides an easy and standard way to do it for module 0 analysis.
(Also, installation requires conda to install. So if you cannot use (or do not want to use) conda for installation, you will want to find a way around this. Or install locally.)
