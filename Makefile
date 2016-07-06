ROOTCONFIG   := root-config
ROOTCINT     := rootcint
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTLIBS     := $(shell $(ROOTCONFIG) --noauxlibs --evelibs)

G4CONFIG     := geant4-config
G4CFLAGS     := $(shell $(G4CONFIG) --cflags) 
G4LDFLAGS    := $(shell $(G4CONFIG) --libs)

GENFITCFLAGS	:=	-I$(GENFIT)/include
GENFITLDFLAGS	:=	-L$(GENFIT)/lib -lgenfit2

CXX           = c++ -g
CXXFLAGS      = -O3 -Wall -fPIC
LD            = c++ -Wall -O3 -g3
LDFLAGS       = -lpthread
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
#LDFLAGS      += $(ROOTLDFLAGS) ${ROOTGLIBS}
LDFLAGS      += $(ROOTLDFLAGS) $(ROOTGLIBS) ${ROOTLIBS}

#CXXFLAGS     += $(G4CFLAGS)
#LDFLAGS      += $(G4LDFLAGS)

#CXXFLAGS     += $(GENFITCFLAGS)
#LDFLAGS      += $(GENFITLDFLAGS)

CXXFLAGS	+=	-I$(OFFLINE_MAIN)/include
LDFLAGS		+=	-L$(OFFLINE_MAIN)/lib -lphg4hit -lg4testbench -lgenfit2exp -lPHGenFit
#LDFLAGS		+=	-L$(OFFLINE_MAIN)/lib -lgenfit2exp -lPHGenFit

CXXFLAGS	+=	-I$(MY_INSTALL)/include
LDFLAGS		+=	-L$(MY_INSTALL)/lib -lgenfit2exp -lPHGenFit

TESTO         = FGEM_PHGenFit_G4Hit_Tracking.o
TEST	        = FGEM_PHGenFit_G4Hit_Tracking

OBJS          = $(TESTO)
PROGRAMS      = $(TEST)

all:            $(PROGRAMS)

.SUFFIXES: .cc .o

$(TEST):       $(TESTO)
	$(LD) $^ -o $@ $(LDFLAGS)
	@echo "$@ done."

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean

clean:
	@echo "Cleanning everything ... "
	@rm $(PROGRAMS) $(OBJS)
