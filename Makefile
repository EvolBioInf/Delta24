SUBDIRS= src docs

.PHONY: subdirs $(SUBDIRS) code-docs clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: src

test: src
	cd test && sh testDeltaComp.sh

code-docs: subdirs

clean:
	cd src && $(MAKE) $@
