#-------------------------------------------------------------------------------
# SuiteSparse_config Makefile
#-------------------------------------------------------------------------------

VERSION = 4.0.2

default: ccode

include SuiteSparse_config.mk

ccode: libsuitesparseconfig.a

all: libsuitesparseconfig.a

library: libsuitesparseconfig.a

libsuitesparseconfig.a: SuiteSparse_config.c SuiteSparse_config.h
	$(CC) $(CF) -c SuiteSparse_config.c
	$(ARCHIVE) libsuitesparseconfig.a SuiteSparse_config.o
	$(RANLIB) libsuitesparseconfig.a
	- $(RM) SuiteSparse_config.o

distclean: purge

purge: clean
	- $(RM) *.o *.a

clean:
	- $(RM) -r $(CLEAN)
