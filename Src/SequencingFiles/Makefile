# Project: ImageExtraction
OUTPUT = SequencingFiles.dll

TOP=..

include $(TOP)/config.mk

SRCS = \
	CIFFile.cs \
	FileHelper.cs \
	IntensityStats.cs \
	Properties/AssemblyInfo.cs \
	RunInfo.cs \
	SerializableToFile.cs
		
DEPS=$(TOP)/Shared/ILMNcommon

all:  deps $(OUTPUT)

install:  deps $(OUTPUT)
	cp $(OUTPUT) $(TOP)/build
	-cp *.mdb $(TOP)/build

deps:
	$(call MAKEDEPS,$(DEPS),install)

REFS = -r:$(TOP)/build/ILMNcommon.dll

$(OUTPUT): $(SRCS)
	gmcs -unsafe $(MCSFLAGS) $^ -target:library -out:$@  $(REFS)
	
clean:
	-rm $(OUTPUT) *.o *.so *.dll *.exe *.mdb; \
	$(call DASH_MAKEDEPS,$(DEPS),clean)