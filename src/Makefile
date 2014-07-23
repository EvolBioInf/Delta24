STATGENDIR= ../lib/libStatGen
LDLIBS= -lStatGen -lz
LDFLAGS= -L $(STATGENDIR)

CXXFLAGS= -I $(STATGENDIR)/include -g -std=c++0x $(LDLIBS) $(LDFLAGS)

Delta24: Delta24.o
	$(CXX) Delta24.o $(CXXFLAGS) -O3 -o Delta24

clean:
	rm -f *.o Delta24
