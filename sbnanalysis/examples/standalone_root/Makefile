CXX=g++
CXXFLAGS=$(shell root-config --cflags) -I$(SBN_ANALYSIS_DIR)/../core
LDFLAGS=$(shell root-config --libs) -L$(SBN_LIB_DIR) -lsbnanalysis_Event

all: eventdump

eventdump: eventdump.cxx
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

