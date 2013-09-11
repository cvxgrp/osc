// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// header file for prox step for multi-period portfolio example

#ifndef PROX_H_GUARD
#define PROX_H_GUARD

#include "../osc.h"
typedef struct {
	double * kappa, * x_term;
}prox_data;

void perturb_data(double x_init[], all_data * data);
void read_in_prox_data(FILE * fb, all_data * data, prox_data * p_data);
void prox(prob_vars * vars, all_data * data, prox_data * p_data);
void free_prox_data(all_data * data,prox_data *p_data);
#endif
