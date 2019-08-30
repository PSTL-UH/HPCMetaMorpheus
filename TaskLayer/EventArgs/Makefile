#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../../Makefile.defs

HEADER =  SingleTaskEventArgs.h     \
	  XmlForTaskListEventArgs.h 


OBJ = 	 SingleTaskEventArgs.o     \
	 XmlForTaskListEventArgs.o 



.SUFFIXES: .o.cpp

all:    
	make $(MMORPHEUS_LIB) 


$(MMORPHEUS_LIB): $(HEADER) $(OBJ) 
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -c $<

clean:  
	rm -rf *.o *~
