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
	cd EngineLayer ; $(MAKE) ; cd ..
	cd TaskLayer ; $(MAKE) ; cd ..
	cd CMD ; $(MAKE) ; cd ..

clean:  
	cd EngineLayer ; $(MAKE) clean ; cd ..
	cd TaskLayer ; $(MAKE) clean ; cd ..
	cd CMD ; $(MAKE) clean ; cd ..
	cd $(MMORPHEUS_LIB_DIR) ; rm -f $(MMORPHEUS_LIB) ; cd ..
	cd $(MMORPHEUS_BIN_DIR) ; rm -f HPCMetaMorpheus ; cd ..	
	rm -rf *.o *~
