SUBDIRS= src

.PHONY:
all: src/Delta24

src/Delta24: src/Delta24.cpp
	cd src && $(MAKE)

test: src/Delta24
	cd test && sh testDeltaComp.sh

clean:
	cd src && $(MAKE) $@
