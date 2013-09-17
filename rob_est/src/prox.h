// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// header file for prox step for robust state estimation example 

#ifndef PROX_H_GUARD
#define PROX_H_GUARD

struct PROX_DATA{
	double M;
};

void prox(prob_vars * vars, all_data * data, prox_data * p_data);
void freeProxData(all_data * data,prox_data *p_data);
#endif
