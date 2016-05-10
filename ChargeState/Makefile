ROOTCONFIG   := root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)


CXX           =
ObjSuf        = o
SrcSuf        = cpp
ExeSuf        =
OutPutOpt     = -o 

CXX           = g++
CXXFLAGS      = -O3 -Wall -fPIC -Wno-deprecated
LD            = g++
LDFLAGS       = -O3

LIBS          = $(SYSLIBS)
GLIBS         = $(SYSLIBS)

LDRFLAGS      = $(LDFLAGS) $(ROOTCFLAGS)

RLIBS         = $(LIBS) $(ROOTGLIBS)

#------------------------------------------------------------------------------
GLOBO       = toi.$(ObjSuf)  chargemodel.$(ObjSuf)
GLOBS       = toi.$(SrcSuf)  chargemodel.$(SrcSuf)

LDFO        = $(GLOBO) ChargeDist.$(ObjSuf)
LDFS        = $(GLOBS) ChargeDist.$(SrcSuf)
LDF         = ChargeDist$(ExeSuf)

DRAWO       = $(GLOBO) chargeDraw.$(ObjSuf)
DRAWS       = $(GLOBS) chargeDraw.$(SrcSuf)
DRAW        = chargeDraw$(ExeSuf)

#OBJS          =  $(LDFO) $(DRAWO)
PROGRAMS      =  $(LDF) $(DRAW)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:           $(PROGRAMS)

$(LDF):        $(LDFO)
	       $(LD) $(LDFLAGS) $(LDFO) $(LIBS) $(OutPutOpt)$@ 
		@echo "$@ DONE"

$(DRAW):       $(DRAWO)
	       $(LD) $(LDRFLAGS) $(DRAWO) $(RLIBS) $(OutPutOpt)$@ 
		@echo "$@ DONE"

chargeDraw.$(ObjSuf): chargeDraw.$(SrcSuf)
	       $(LD) $(CXXFLAGS) $(LDRFLAGS) -c chargeDraw.$(SrcSuf) $(OutPutOpt)$@ 

clean:
		@rm -f $(LDFO) $(DRAWO) core

distclean:      clean
		@rm -f $(PROGRAMS) $(LDFO)

