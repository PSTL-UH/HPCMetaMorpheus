#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../../Makefile.defs

HEADER =  AllowedIntervalWithNotch.h                 \
	  IntervalMassDiffAcceptor.h                 \
	  MassDiffAcceptor.h                         \
	  SingleAbsoluteAroundZeroMassDiffAcceptor.h \
	  DotMassDiffAcceptor.h                      \
	  LowTheoreticalDiffAcceptor.h               \
	  OpenMassDiffAcceptor.h                     \
	  SinglePpmAroundZeroMassDiffAcceptor.h


OBJ = 	  AllowedIntervalWithNotch.o                 \
	  MassDiffAcceptor.o                         \
	  DotMassDiffAcceptor.o                      \
	  IntervalMassDiffAcceptor.o                 \
	  SingleAbsoluteAroundZeroMassDiffAcceptor.o \
	  LowTheoreticalDiffAcceptor.o               \
	  OpenMassDiffAcceptor.o                     \
	  SinglePpmAroundZeroMassDiffAcceptor.o


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
