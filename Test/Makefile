#
# Copyright (c) 2019-2020   University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../Makefile.defs

CEXECS = BinGenerationTest \
	 TestToml




.SUFFIXES: .cxx .h  

###################################################################################################

all: cxxtests 

cxxtests: $(CEXECS)


###################################################################################################
### CXXTESTS                                                                                      ###
###################################################################################################
TestToml: TestToml.o 
	$(CXX) $^ -o $@ $(LDFLAGS)

BinGenerationTest: BinGenerationTest.o 
	$(CXX) $(CXXFLAGS) -c TestDataFile.cpp
	$(CXX) $^ -o $@ TestDataFile.o $(LDFLAGS)

###################################################################################################


clean:
	rm -rf *.o *~

	rm -f TestToml BinGenerationTest
