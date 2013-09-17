// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// header file for both osc.out warm_start.out
// this file defines structures that will be used throughout

#ifndef OSC_H_GUARD                                                              
#define OSC_H_GUARD

#include <stdbool.h>

// struct that containing standard problem data
typedef struct PROBLEM_DATA {
	size_t n,m,T,nc;
	int LNZ;
	double rho;
	double * x_init;
	int * Lp, *Lr; // factorization of KKT matrix
	double * Lx, *D; // factorization of KKT matrix
	double * RHS; //component of RHS of KKT sys that never changes
	int *P,*Pinv; // permutation of KKT matrix for factorization
	double alpha; //relaxation parameter;
} all_data;

// struct containing problem variables
// we are consistent with the notation in the paper
// x_t and u_t correspond to \tilde x and \tilde u resp.
typedef struct PROBLEM_VARS{
	double *x,*u,*x_t,*u_t,*z,*y;
} prob_vars;

// structure used to report timing
typedef struct TT{
	double prox_time, lin_sys_time, total_time, itns;
} Timings;

typedef struct PROX_DATA prox_data;

// function prototypes:
// OSC:
Timings osc(prob_vars * vars, all_data * data, prox_data * p_data);

// test convergence of residuals
bool testConvergence(all_data * data,prob_vars * vars,double *x_t_old,double *u_t_old,double*rn,double*dn);

// solves the KKT system
void solveLinSys(prob_vars * vars,all_data *data);

// update the (scaled) dual variables 
void updateDualVars(prob_vars * vars,all_data *data);

// simple functions to add or subtract two arrays and store result in a
void addArray(double * a,double *b,int len);
void subArray(double * a,double *b,int len);

// helper function to print all vars for debugging
void printAll(all_data * data, prob_vars * vars);

// performs the relaxation of the iterates as described in the paper
void relax(prob_vars * vars, all_data* data);

// frees memory associated with data and variables at end of program
void freeData(all_data* data);
void freeVars(prob_vars * vars);

void prox(prob_vars * vars, all_data * data, prox_data * p_data);
#endif
