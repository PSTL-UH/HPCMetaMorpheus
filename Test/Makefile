#
# Copyright (c) 2019-2020   University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../Makefile.defs

CEXECS = BinGenerationTest  \
	 CalibrationTests   \
	 XLTest             \
	 XLSearchOutputTest \
	 SearchEngineTests  \
	 IndexEngineTests   \
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

CalibrationTests: CalibrationTests.o 
	$(CXX) $^ -o $@ $(LDFLAGS)

XLTest: XLTest.o 
	$(CXX) $^ -o $@ $(LDFLAGS)

XLSearchOutputTest: XLSearchOutputTest.o 
	$(CXX) $^ -o $@ $(LDFLAGS)

SearchEngineTests: SearchEngineTests.o 
	$(CXX) $^ -o $@ $(LDFLAGS)

IndexEngineTests: IndexEngineTests.o 
	$(CXX) $^ -o $@ $(LDFLAGS)

BinGenerationTest: BinGenerationTest.o ../lib/libmmorpheus.a
	$(CXX) $(CXXFLAGS) -c TestDataFile.cpp
	$(CXX) $^ -o $@ TestDataFile.o $(LDFLAGS)

###################################################################################################


clean:
	rm -rf *.o *~

	rm -f TestToml BinGenerationTest CalibrationTests XLTest XLSearchOutputTest
	rm -f SearchEngineTests IndexEngineTest
	rm -f TestProteinSplitAcrossFiles.xml TestProteinSplitAcrossFiles1.ind
	rm -f TestProteinSplitAcrossFiles1.mzML TestProteinSplitAcrossFiles2.mzML
	rm -f pep.XML.pep.XM BinGenerationTest.xml BinGenerationTest.mzML
