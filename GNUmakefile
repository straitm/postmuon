LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_novadaq.mk
include SoftRelTools/arch_spec_novadb.mk

override CPPFLAGS := -I$(NOVADAQ_INC) $(CPPFLAGS) -I$(NUTOOLS_INC) 

LIBLINK     := \
-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) \
-l$(PACKAGE)

override LIBLIBS += \
$(LOADLIBES) \
-L$(ART_LIB) \
-L$(NOVADAQ_LIB) \
-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) \
-lDAQDataFormats -lCalibrator -lRawData
#-lChannelInfo -lRecoBase 
