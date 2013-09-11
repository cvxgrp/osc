// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// header file for prox step for supply chain example

#ifndef PROX_H_GUARD
#define PROX_H_GUARD

#include "../osc.h"

// structure declarations	
typedef struct {
	unsigned int numsource;
	unsigned int * idx_source;
	unsigned int ** connections;
	double C;
	double U;
}prox_data;

// function prototypes
void perturb_data(double x_init[], all_data * data);
void read_in_prox_data(FILE * fb, all_data * data, prox_data * p_data);
void prox(prob_vars * vars, all_data * data, prox_data * p_data);
void free_prox_data(all_data * data, prox_data *p_data);
#endif

