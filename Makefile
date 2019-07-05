#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include Makefile.defs


all:    
	cd EngineLayer ; make ; cd ..
	cd TaskLayer ; make ; cd ..
	cd CMD ; make ; cd ..

clean:  
	cd EngineLayer ; make clean ; cd ..
	cd TaskLayer ; make clean ; cd ..
	cd CMD ; make clean ; cd ..
	cd $(MMORPHEUS_LIB_DIR) ; rm -f $(MMORPHEUS_LIB) ; cd ..
	rm -rf *.o *~
