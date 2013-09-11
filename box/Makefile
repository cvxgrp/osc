# A simple makefile for managing build of project composed of C source files.
#
UNAME = $(shell uname -s)
# It is likely that default C compiler is already gcc, but explicitly
# set, just to be sure
CC = gcc

# The CFLAGS variable sets compile flags for gcc:
#  -g          compile with debug information
#  -Wall       give all diagnostic warnings
#  -pedantic   require compliance with ANSI standard
#  -O0         do not optimize generated code
#  -std=gnu99  use the Gnu C99 standard language definition
#  -m32        emit code for IA32 architecture

CFLAGS = -g -Wall -pedantic -O3 -fopenmp 

# The LDFLAGS variable sets flags for linker
LDFLAGS = -lm -lldl -lamd -L ../AMD/Lib/ -L ../LDL/Lib/

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

# install locations for packages: 
CFLAGS += -I ../AMD/Include -I ../LDL/Include -I ../UFconfig/

# In this section, you list the files that are part of the project.
# If you add/change names of header/source files, here is where you
# edit the Makefile.

HEADERS = ../cholesky.h ../osc.h prox.h
SOURCES = ../cholesky.c prox.c osc.c warm_start.c
OBJECTS = $(SOURCES:.c=.o)
TARGETS = osc.out warm_start.out

# used when make is invoked with no argument. Given the definitions
# above, this Makefile file will build the one named TARGET and
# assume that it depends on all the named OBJECTS files.

default: $(TARGETS)

osc.out : ../cholesky.o osc.o prox.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

warm_start.out : ../cholesky.o warm_start.o prox.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# In make's default rules, a .o automatically depends on its .c file
# (so editing the .c will cause recompilation into its .o file).
# The line below creates additional dependencies, most notably that it
# will cause the .c to reocmpiled if any included .h file changes.

Makefile.dependencies:: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -MM $(SOURCES) > Makefile.dependencies

-include Makefile.dependencies

# Phony means not a "real" target, it doesn't build anything
# The phony target "clean" that is used to remove all compiled object files.

.PHONY: clean

clean:
	@rm -rf $(TARGETS) $(OBJECTS) core Makefile.dependencies *.o *.a
