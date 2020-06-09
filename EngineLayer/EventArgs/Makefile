#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../../Makefile.defs

HEADER =  MyRecursiveEventArgs.h          \
	  SingleEngineEventArgs.h         \
	  SingleFileEventArgs.h           \
	  StringListEventArgs.h           \
	  ProgressEventArgs.h             \
	  SingleEngineFinishedEventArgs.h \
	  StringEventArgs.h


OBJ = 	  MyRecursiveEventArgs.o          \
	  SingleEngineEventArgs.o         \
	  SingleFileEventArgs.o           \
	  StringListEventArgs.o           \
	  ProgressEventArgs.o             \
	  SingleEngineFinishedEventArgs.o \
	  StringEventArgs.o



.SUFFIXES: .o.cpp

all:    
	$(MAKE) $(MMORPHEUS_LIB) 


$(MMORPHEUS_LIB): $(HEADER) $(OBJ) 
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -c $<

clean:  
	rm -rf *.o *~
