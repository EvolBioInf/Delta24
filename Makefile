SUBDIRS= src docs

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: src

test: src
	cd test && sh testDeltaComp.sh

code-docs: subdirs

.PHONY: clean
clean:
	cd src && $(MAKE) $@
