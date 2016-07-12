HDF5FLAGS    := -I/cvmfs/oasis.opensciencegrid.org/osg/modules/hdf5/1.8.9/include
HDF5LIBS     := -L/cvmfs/oasis.opensciencegrid.org/osg/modules/hdf5/1.8.9/lib

ROOTCFLAGS   := -pthread -m64 -I/cvmfs/oasis.opensciencegrid.org/osg/modules/root/5.34-21/include/root
ROOTLIBS     := -L/cvmfs/oasis.opensciencegrid.org/osg/modules/root/5.34-21/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -Wl,-rpath=/cvmfs/oasis.opensciencegrid.org/osg/modules/hdf5/1.8.9/lib/
LD            = g++
LDFLAGS       = -O3

CXXFLAGS     += $(ROOTCFLAGS) $(HDF5FLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) $(HDF5LIBS) -lhdf5 -lhdf5_cpp

vpath %.cpp src
objdir     = obj

SRC        = cll.cpp eos.cpp eos_hotqcd.cpp trancoeff.cpp fld.cpp hdo.cpp main.cpp mcgstream.cpp rmn.cpp cornelius.cpp 
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = hlle_visc
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
	
obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
