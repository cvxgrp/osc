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

// function prototypes:

// simple function to open supplied file
int open_file(int argc, char ** argv, int idx, char * default_file, FILE ** fb);

// read in the data from supplied data file into all_data struct
int read_in_data(FILE * fp,all_data * data);

// initialize the variables to zero (with correct memory allocation)
void init_vars(prob_vars * vars,all_data* data);

// test convergence of residuals
bool test_convergence(all_data * data,prob_vars * vars,double *x_t_old,double *u_t_old,double*rn,double*dn);

// solves the KKT system
void solve_lin_sys(prob_vars * vars,all_data *data);

// update the (scaled) dual variables 
void update_dual_vars(prob_vars * vars,all_data *data);

// simple functions to add or subtract two arrays and store result in a
void add_array(double * a,double *b,int len);
void sub_array(double * a,double *b,int len);

// helper function to print all vars for debugging
void print_all(all_data * data, prob_vars * vars);

// performs the relaxation of the iterates as described in the paper
void relax(prob_vars * vars, all_data* data);

// frees memory associated with data and variables at end of program
void free_data(all_data* data);
void free_vars(prob_vars * vars);

#endif
