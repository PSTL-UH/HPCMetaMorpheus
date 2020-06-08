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
	make $(MMORPHEUS_LIB) 
	cd CalibrationTask ; make ; cd ..
	cd EventArgs ; make ; cd ..
	cd SearchTask ; make ; cd ..
	cd XLSearchTask ; make ; cd ..

#	cd GPTMDTask ; make ; cd ..


$(MMORPHEUS_LIB): $(HEADER) $(OBJ) 
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -c $<

clean:  
	rm -rf *.o *~
	cd CalibrationTask ; make clean ; cd ..
	cd EventArgs ; make clean ; cd ..
	cd SearchTask ; make clean ; cd ..
	cd XLSearchTask ; make clean ; cd ..

#	cd GPTMDTask ; make clean ; cd ..
