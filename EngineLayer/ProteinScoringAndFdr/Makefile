#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../../Makefile.defs

HEADER =  FdrCategory.h                \
	  FdrClassifier.h              \
	  ProteinScoringAndFdrEngine.h \
	  ProteinScoringAndFdrResults.h


OBJ = 	  FdrCategory.o                \
	  FdrClassifier.o              \
	  ProteinScoringAndFdrEngine.o \
	  ProteinScoringAndFdrResults.o



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
