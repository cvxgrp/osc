# A Splitting Method for Optimal Control
# by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
# ===========================================================

UNAME = $(shell uname -s)
CC = gcc

CFLAGS = -g -Wall -pedantic -O3 -fopenmp 

LDFLAGS = -lm src/external/LDL/Lib/libldl.a src/external/AMD/Lib/libamd.a 

HEADERS = src/osc.h src/cholesky.h src/run_osc.h
SOURCES = $(HEADERS:.h=.c)
OBJECTS = $(HEADERS:.h=.o) 
TARGETS = box/run_osc finance/run_osc rob_est/run_osc sup_ch/run_osc

ifeq ($(UNAME), Darwin)
	CFLAGS   += -std=c99
else
	CFLAGS   += -std=gnu99
endif

CFLAGS += -Isrc/external/AMD/Include -Isrc/external/LDL/Include -Isrc/external/SuiteSparse_config

default: packages examples

packages :
	cd src/external/AMD && make
	cd src/external/LDL && make

examples : box finance rob_est sup_ch

box : $(OBJECTS) box/src/prox.o
	$(CC) $(CFLAGS) -o box/run_osc $^ $(LDFLAGS)

finance : $(OBJECTS) finance/src/prox.o
	$(CC) $(CFLAGS) -o finance/run_osc $^ $(LDFLAGS)

rob_est : $(OBJECTS) rob_est/src/prox.o
	$(CC) $(CFLAGS) -o rob_est/run_osc $^ $(LDFLAGS)

sup_ch : $(OBJECTS) sup_ch/src/prox.o
	$(CC) $(CFLAGS) -o sup_ch/run_osc $^ $(LDFLAGS)

.PHONY: clean

clean:
	@rm -rf $(TARGETS) $(OBJECTS) core Makefile.dependencies *.o *.a box/*.o finance/*.o rob_est/*.o sup_ch/*.o
	cd src/external/AMD && make purge
	cd src/external/LDL && make purge
