#ifndef OSC_H_GUARD                                                              
#define OSC_H_GUARD

#include <stdbool.h>
typedef struct PROBLEM_DATA {
	size_t n,m,T,nc;
	int LNZ;
	double rho;
	double * x_init;
	int * Lp, *Lr;
	double * Lx, *D;
	double * RHS; //component of RHS of KKT sys that never changes
	int *P,*Pinv;
	double alpha; //relaxation parameter;
} all_data;

typedef struct PROBLEM_VARS{
	double *x,*u,*x_t,*u_t,*z,*y;
} prob_vars;

int read_in_data(FILE * fp,all_data * data);
void update_dual(prob_vars * vars,all_data * data);
void init_vars(prob_vars * vars,all_data* data);
bool test_convergence(all_data * data,prob_vars * vars,double *x_t_old,double *u_t_old,double*rn,double*dn);
void solve_lin_sys(prob_vars * vars,all_data *data);
void update_dual_vars(prob_vars * vars,all_data *data);
void add_array(double * a,double *b,int len);
void sub_array(double * a,double *b,int len);
#endif
