LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lart_Framework_Services_Registry -lPhysics -lart_Persistency_Common -lart_Framework_Core -l$(PACKAGE)

override CPPFLAGS += -I$(NOVARWGT_INC)


include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_ifdhart.mk
include SoftRelTools/arch_spec_ifdhc.mk
include SoftRelTools/arch_spec_nutools.mk
include SoftRelTools/arch_spec_genie.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB)  -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) \
-lRecoBase \
-lStandardRecord \
-lStandardRecordProxy \
-lCAFAna -lCAFAnaVars -lCAFAnaCuts -lCAFAnaCore \
 -lMCCheater -lFuzzyKVertex \
-lMinuit2 \
-lRecoJMShower \
-lShowerLID \
-lBackTracker_service \
-lMCReweight
