# SVJ
Pythia8.309 is used below, and compiled using [Makefile](Makefile).<br />
[sift.h](sift.h) is the header file which contains the code to cluster the event, and gives list of jets as output.
[QCD](QCD.dat) and [HV](HV.dat) are the pythia command files.
[tutorial1.cc](tutorial1.cc) is a sample pythia program which generates events, and then forms jets and subjets.
Caurrently I am trying to make the program faster by using VP Trees to find the closest pair of particles.
