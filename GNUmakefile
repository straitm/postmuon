LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_root.mk

override CPPFLAGS := $(CPPFLAGS) -I$(NOVADAQ_INC)

#LIBLINK     := \
#-l$(PACKAGE)
#-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
#-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) \

#override LIBLIBS += \
#$(LOADLIBES) \
#-L$(ART_LIB) \
#-L$(NOVADAQ_LIB) \
#-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
#-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR)
