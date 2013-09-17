// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// header file for prox step for supply chain example

#ifndef PROX_H_GUARD
#define PROX_H_GUARD

// structure declarations	
struct PROX_DATA{
	unsigned int numsource;
	unsigned int * idx_source;
	unsigned int ** connections;
	double C;
	double U;
};

// function prototypes
void prox(prob_vars * vars, all_data * data, prox_data * p_data);
void freeProxData(all_data * data, prox_data *p_data);
#endif

