STATGENDIR= ../lib/libStatGen
LDLIBS= -lStatGen -lz
LDFLAGS= -L $(STATGENDIR)

CXXFLAGS= -I $(STATGENDIR)/include -g -std=c++11 -Wall -W $(LDLIBS) $(LDFLAGS) -O3 -fstack-protector -fopenmp

Delta24: Delta24.o mapNucl.o dml.o nucl.o
	$(CXX) $^ $(CXXFLAGS) -o Delta24

.PHONY: prof
prof: CXXFLAGS+=-pg
prof: clean Delta24
	./Delta24 ../test/test2.bam 1 6 1
	gprof Delta24 > profile

.PHONY: debug
debug: CXXFLAGS= -I $(STATGENDIR)/include -g -std=c++11 -Wall -W $(LDLIBS) $(LDFLAGS) -Og -fstack-protector-all
debug: Delta24

clean:
	rm -f *.o Delta24
