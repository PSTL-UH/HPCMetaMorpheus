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

LINKLIBS = $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB)

.SUFFIXES: .o.cpp

all:
	make HPCMetaMorpheus

$(OBJ):	
.cpp.o:	
	$(CXX) $(CXXFLAGS) -c $<

HPCMetaMorpheus: $(OBJ) $(LINKLIBS)
	$(CXX) $(CXXFLAGS) $^ -o $(MMORPHEUS_BIN_DIR)/$(MMORPHEUS_EXEC) $(LDFLAGS)


clean:  
	rm -rf *.o *~ HPCMetaMorpheus
