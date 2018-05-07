LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_ifdhart.mk
include SoftRelTools/arch_spec_ifdhc.mk
include SoftRelTools/arch_spec_nutools.mk
include SoftRelTools/arch_spec_genie.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB)
