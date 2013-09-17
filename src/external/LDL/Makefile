#------------------------------------------------------------------------------
# LDL Makefile
#------------------------------------------------------------------------------

VERSION = 2.1.0

default: library

include ../SuiteSparse_config/SuiteSparse_config.mk

# compile just the C-callable libraries (not Demos)
library:
	( cd Lib    ; $(MAKE) )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Lib    ; $(MAKE) clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd Lib    ; $(MAKE) purge )
