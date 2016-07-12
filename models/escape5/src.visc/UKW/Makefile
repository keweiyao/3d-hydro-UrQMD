CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC $(ROOTCFLAGS) -g
LD            = g++
LDFLAGS       = $(ROOTLIBS) -lgfortran -g
SOFLAGS       = -shared
FFLAGS        = -fPIC $(ROOTCFLAGS) -g


LIBS          = $(SYSLIBS)
GLIBS         = $(SYSLIBS)

OBJCXX      = main.o
OBJS        = mainf.o a.o siglookup.o addpart.o init.o output.o string.o input.o tabinit.o paulibl.o angdis.o dectim.o \
     proppot.o ukw.o anndec.o delpart.o getmass.o upmerge.o detbal.o getspin.o iso.o upyth.o blockres.o ityp2pdg.o \
     dwidth.o jdecay2.o whichres.o boxprg.o hepchg.o cascinit.o hepcmp.o make22.o hepnam.o  \
     coload.o numrec.o saveinfo.o scatter.o error.o
              
TARGET	    = test
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS) $(OBJCXX)
		$(LD) $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(OBJCXX) $(TARGET)

$(OBJCXX) : %.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $<
	
$(OBJS) : %.o : %.f
	$(F) $(FFLAGS) -c $<
