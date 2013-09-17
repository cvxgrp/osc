# Produced by CVXGEN, 2012-08-15 14:12:46 -0700.
# CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com.
# The code in this file is Copyright (C) 2006-2011 Jacob Mattingley.
# CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
# applications without prior written permission from Jacob Mattingley.

# Filename: Makefile.
# Description: Basic Makefile.

OPT = -Wall

# for valgrind.
#OPT=-g -O0

OPT=-Os

# libmath is needed for sqrt. can remove in final product.
LDLIBS = -lm -lldl -L LDL/Lib/
CFLAGS = $(OPT) -I LDL/Include -I UFconfig/
CC = gcc

.PHONY: all
all: testsolver

packages:
	cd UFconfig && make purge && make
	cd LDL && make purge && make

testsolver: solver.o matrix_support.o ldl.o testsolver.o util.o cholesky.o

# Include util.o for random functions and easy matrix printing.
#testsolver: solver.o matrix_support.o ldl.o util.o testsolver.o

solver.o: solver.h
matrix_support.o: solver.h
	$(CC) -O0 -c -o matrix_support.o matrix_support.c
ldl.o: solver.h cholesky.h
	$(CC) -O0 -I LDL/Include -I UFconfig/ -c -o ldl.o ldl.c
util.o: solver.h
testsolver.o: solver.h

.PHONY : clean
clean :
	-rm -f *.o testsolver
	cd UFconfig && make purge
	cd LDL && make purge
