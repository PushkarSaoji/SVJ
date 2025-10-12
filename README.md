# SVJ
Pythia8.309 is used. The `Makefile.inc` from compiling `Pythia8` should be copied here.
The test code is compiled using [Makefile](Makefile).<br />
[sift.h](sift.h) is the header file that contains the code to cluster the event as per the SIFT(Scale-Invariant Filter Tree) algorithm, and gives a list of jets as output.
[QCD](QCD.dat) and [HV](HV.dat) are the Pythia command files to generate QCD and HV (Hidden Valley) events, respectively.
[tutorial1.cc](tutorial1.cc) is a sample Pythia program that generates events and then forms jets and subjets.
In [tutorial2.cc](tutorial2.cc), I am trying to make the program faster by using VP Trees to find the closest pair of particles.
[tutorial3.cc](tutorial3.cc) contains the implementation of the Heap-Based algorithm, achieving significant speedup over the other implementations. 
In [tutorial4.cc](tutorial4.cc), ghost particles are added to check the time complexity for larger events.
Data for time taken vs number of particles can be found at [t_c.csv](t_c.csv).
