#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../Makefile.defs

HEADER = DbForTask.h              \
	 EverythingRunnerEngine.h \
	 FileSpecificParameters.h \
	 MetaMorpheusTask.h       \
	 MyFileManager.h          \
	 MyTaskResults.h          \
	 PepXMLWriter.h


OBJ = 	 DbForTask.o              \
	 PepXMLWriter.o           \
	 MyFileManager.o          \
	 MyTaskResults.o          \
	 FileSpecificParameters.o \
	 MetaMorpheusTask.o       \
	 EverythingRunnerEngine.o 


.SUFFIXES: .o.cpp

all:    
	$(MAKE) $(MMORPHEUS_LIB) 
	cd CalibrationTask ; $(MAKE) ; cd ..
	cd EventArgs ; $(MAKE) ; cd ..
	cd SearchTask ; $(MAKE) ; cd ..
	cd XLSearchTask ; $(MAKE) ; cd ..

#	cd GPTMDTask ; $(MAKE) ; cd ..


$(MMORPHEUS_LIB): $(HEADER) $(OBJ) 
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -c $<

clean:  
	rm -rf *.o *~
	cd CalibrationTask ; $(MAKE) clean ; cd ..
	cd EventArgs ; $(MAKE) clean ; cd ..
	cd SearchTask ; $(MAKE) clean ; cd ..
	cd XLSearchTask ; $(MAKE) clean ; cd ..

#	cd GPTMDTask ; make clean ; cd ..
