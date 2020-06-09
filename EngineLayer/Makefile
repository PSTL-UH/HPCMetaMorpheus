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
	$(MAKE) $(MMORPHEUS_LIB) 
	cd EventArgs; $(MAKE) ; cd ..
	cd Indexing; $(MAKE) ; cd ..
	cd FdrAnalysis; $(MAKE) ; cd ..
	cd ModernSearch; $(MAKE) ; cd ..
	cd CrosslinkSearch; $(MAKE) ; cd ..
	cd PrecursorSearchModes; $(MAKE) ; cd ..
	cd Calibration; $(MAKE) ; cd ..
	cd ClassicSearch; $(MAKE) ; cd ..
	cd HistogramAnalysis; $(MAKE) ; cd ..
	cd Localization; $(MAKE) ; cd ..
	cd ModificationAnalysis; $(MAKE) ; cd ..
	cd NonSpecificEnzymeSearch; $(MAKE) ; cd ..
	cd ProteinParsimony; $(MAKE) ; cd ..
	cd ProteinScoringAndFdr; $(MAKE) ; cd ..


#	cd Gptmd; $(MAKE) ; cd ..
#	cd MetaDraw; $(MAKE) ; cd ..


$(MMORPHEUS_LIB): $(HEADER) $(OBJ) 
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -c $<

clean:  
	rm -rf *.o *~
	cd EventArgs; $(MAKE) clean; cd ..
	cd Indexing; $(MAKE) clean; cd ..
	cd FdrAnalysis; $(MAKE) clean; cd ..
	cd ModernSearch; $(MAKE) clean; cd ..
	cd CrosslinkSearch; $(MAKE) clean; cd ..
	cd PrecursorSearchModes; $(MAKE) clean; cd ..
	cd Calibration; $(MAKE) clean; cd ..
	cd ClassicSearch; $(MAKE) clean; cd ..
	cd HistogramAnalysis; $(MAKE) clean; cd ..
	cd Localization; $(MAKE) clean; cd ..
	cd ModificationAnalysis; $(MAKE) clean; cd ..
	cd NonSpecificEnzymeSearch; $(MAKE) clean; cd ..
	cd ProteinParsimony; $(MAKE) clean; cd ..
	cd ProteinScoringAndFdr; $(MAKE) clean; cd ..

#	cd Gptmd; $(MAKE) clean; cd ..
#	cd MetaDraw; $(MAKE) clean; cd ..
