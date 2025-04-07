# Include the configuration.
-include Makefile.inc

CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl
PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)

#tutorial%: tutorial%.cc
#	g++ -I/home/pushkar/Downloads/pythia8309/include $@.cc -o $@ -lpythia8 -L/home/pushkar/Downloads/pythia8307/lib
	# Examples without external dependencies.
tutorial%: $(PYTHIA) tutorial%.cc
	$(CXX) $@.cc -o $@ $(CXX_COMMON) $(HEPMC3_LIB) -Wl,-E -Wl,--hash-style=gnu -fno-inline -fno-omit-frame-pointer

