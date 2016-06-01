#
# relies on IDS definitions
#

IDS := PYFLOW_ID

PYFLOW_ID := pyflow-1.1.8

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

id-to-dir = $($1)
DIRS := $(foreach I,$(IDS), $(call id-to-dir,$I))

TARGET_LINKS := pyflow

.PHONY: $(IDS)
all: $(TARGET_LINKS) 

$(IDS): $(DIRS)

safe-expand = \
  mkdir -p $1.extract.tmp; \
  $2 -dc < $1 | tar -xf - -C $1.extract.tmp; \
  mv $1.extract.tmp/* .; \
  rmdir $1.extract.tmp


maybe-expand = $(if $(wildcard $1),$(call safe-expand,$1,$2,$3),)

$(DIRS):
	$(call maybe-expand,$@.tar.bz2,bzip2,$@)
	$(call maybe-expand,$@.tbz,bzip2,$@)
	$(call maybe-expand,$@.tar.gz,gzip,$@)
	$(call maybe-expand,$@.tgz,gzip,$@)

pyflow: $(PYFLOW_ID)
	ln -sf $< $@

clean:
	$(RM) -rf $(DIRS) $(TARGET_LINKS) 

