#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../Makefile.defs

HEADER = CommonParameters.h           \
	 GlobalSettings.h             \
	 GlobalVariables.h            \
	 IGlobalSettings.h            \
	 IScan.h                      \
	 MetaMorpheusEngine.h         \
	 MetaMorpheusEngineResults.h  \
	 MetaMorpheusException.h      \
	 Ms2ScanWithSpecificMass.h    \
	 PeptideSpectralMatch.h       \
	 ScanWithIndexAndNotchInfo.h



OBJ = 	 CommonParameters.o           \
	 GlobalSettings.o             \
	 GlobalVariables.o            \
	 IGlobalSettings.o            \
	 IScan.o                      \
	 MetaMorpheusEngine.o         \
	 MetaMorpheusEngineResults.o  \
	 MetaMorpheusException.o      \
	 Ms2ScanWithSpecificMass.o    \
	 PeptideSpectralMatch.o       \
	 ScanWithIndexAndNotchInfo.o


.SUFFIXES: .o.cpp

all:    
	make $(MMORPHEUS_LIB) 
	cd Calibration; make ; cd ..
	cd ClassicSearch; make ; cd ..
	cd CrosslinkSearch; make ; cd ..
	cd EventArgs; make ; cd ..
	cd FdrAnalysis; make ; cd ..
	cd Gptmd; make ; cd ..
	cd HistogramAnalysis; make ; cd ..
	cd Indexing; make ; cd ..
	cd Localization; make ; cd ..
	cd MetaDraw; make ; cd ..
	cd ModernSearch; make ; cd ..
	cd ModificationAnalysis; make ; cd ..
	cd NonSpecificEnzymeSearch; make ; cd ..
	cd PrecursorSearchModes; make ; cd ..
	cd ProteinParsimony; make ; cd ..
	cd ProteinScoringAndFdr; make ; cd ..


$(MMORPHEUS_LIB): $(HEADER) $(OBJ) 
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -c $<

clean:  
	rm -rf *.o *~
	cd Calibration; make clean; cd ..
	cd ClassicSearch; make clean; cd ..
	cd CrosslinkSearch; make clean; cd ..
	cd EventArgs; make clean; cd ..
	cd FdrAnalysis; make clean; cd ..
	cd Gptmd; make clean; cd ..
	cd HistogramAnalysis; make clean; cd ..
	cd Indexing; make clean; cd ..
	cd Localization; make clean; cd ..
	cd MetaDraw; make clean; cd ..
	cd ModernSearch; make clean; cd ..
	cd ModificationAnalysis; make clean; cd ..
	cd NonSpecificEnzymeSearch; make clean; cd ..
	cd PrecursorSearchModes; make clean; cd ..
	cd ProteinParsimony; make clean; cd ..
	cd ProteinScoringAndFdr; make clean; cd ..
