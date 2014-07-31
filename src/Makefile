STATGENDIR= ../lib/libStatGen
LDLIBS= -lStatGen -lz
LDFLAGS= -L $(STATGENDIR)

CXXFLAGS= -I $(STATGENDIR)/include -g -std=c++0x -Wall -W $(LDLIBS) $(LDFLAGS) -O3

Delta24: Delta24.o matHash.o bam24.o
	$(CXX) $^ $(CXXFLAGS) -o Delta24

clean:
	rm -f *.o Delta24
