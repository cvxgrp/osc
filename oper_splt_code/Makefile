# A Splitting Method for Optimal Control
# by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
# ===========================================================

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
#-m32
#-std=gnu99

# The LDFLAGS variable sets flags for linker
LDFLAGS = -lm -lldl -lamd -L AMD/Lib/ -L LDL/Lib/
#-lgsl -lgslcblas -m32

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

# install locations for packages: 
CFLAGS += -I AMD/Include -I LDL/Include -I UFconfig/

# In this section, you list the files that are part of the project.
# If you add/change names of header/source files, here is where you
# edit the Makefile.

# used when make is invoked with no argument. Given the definitions
# above, this Makefile file will build the one named TARGET and
# assume that it depends on all the named OBJECTS files.

default: clean copy_c_files packages examples

copy_c_files :
	cp osc.c box
	cp osc.c finance
	cp osc.c rob_est
	cp osc.c sup_ch
	cp warm_start.c box
	cp warm_start.c finance
	cp warm_start.c rob_est
	cp warm_start.c sup_ch

packages :
	cd UFconfig && make
	cd AMD && make
	cd LDL && make

examples :
	cd box && make
	cd finance && make
	cd rob_est && make
	cd sup_ch && make

# In make's default rules, a .o automatically depends on its .c file
# (so editing the .c will cause recompilation into its .o file).
# The line below creates additional dependencies, most notably that it
# will cause the .c to reocmpiled if any included .h file changes.

# Phony means not a "real" target, it doesn't build anything
# The phony target "clean" that is used to remove all compiled object files.

.PHONY: clean

clean:
	@rm -rf $(TARGETS) $(OBJECTS) core Makefile.dependencies *.o *.a
	cd AMD && make purge
	cd LDL && make purge
	cd UFconfig && make purge
	cd box && make clean
	cd finance && make clean
	cd rob_est && make clean
	cd sup_ch && make clean
