#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../Makefile.defs

OBJ = 	 Program.o

.SUFFIXES: .o.cpp


$(OBJ):	
.cpp.o:	
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(CXXFLAGS) $^ -o $(MMORPHEUS_BIN_DIR)/HPCMetaMorpheus $(LDFLAGS)

clean:  
	rm -rf *.o *~ HPCMetaMorpheus
