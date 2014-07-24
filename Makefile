SUBDIRS= src

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: src

test: src
	cd test && sh testDeltaComp.sh

clean: src
