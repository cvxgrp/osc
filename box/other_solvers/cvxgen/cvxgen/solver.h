/* Produced by CVXGEN, 2012-08-15 14:12:43 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.h. */
/* Description: Header file with relevant definitions. */

#ifndef SOLVER_H
#define SOLVER_H

/* Uncomment the next line to remove all library dependencies. */
/*#define ZERO_LIBRARY_MODE */

#ifdef MATLAB_MEX_FILE
/* Matlab functions. MATLAB_MEX_FILE will be defined by the mex compiler. */
/* If you are not using the mex compiler, this functionality will not intrude, */
/* as it will be completely disabled at compile-time. */
#include "mex.h"
#else
#ifndef ZERO_LIBRARY_MODE
#include <stdio.h>
#endif
#endif

/* Space must be allocated somewhere (testsolver.c, csolve.c or your own */
/* program) for the global variables vars, params, work and settings. */
/* At the bottom of this file, they are externed. */

#ifndef ZERO_LIBRARY_MODE
#include <math.h>
#define pm(A, m, n) printmatrix(#A, A, m, n, 1)
#endif

typedef struct Params_t {
  double x_0[50];
  double Q[2500];
  double R[400];
  double A[2500];
  double B[1000];
  double u_max[1];
  double u_min[1];

  double *x[1];
} Params;

typedef struct Vars_t {
  double *u_0; /* 20 rows. */
  double *x_1; /* 50 rows. */
  double *u_1; /* 20 rows. */
  double *x_2; /* 50 rows. */
  double *u_2; /* 20 rows. */
  double *x_3; /* 50 rows. */
  double *u_3; /* 20 rows. */
  double *x_4; /* 50 rows. */
  double *u_4; /* 20 rows. */
  double *x_5; /* 50 rows. */
  double *u_5; /* 20 rows. */
  double *x_6; /* 50 rows. */
  double *u_6; /* 20 rows. */
  double *x_7; /* 50 rows. */
  double *u_7; /* 20 rows. */
  double *x_8; /* 50 rows. */
  double *u_8; /* 20 rows. */
  double *x_9; /* 50 rows. */
  double *u_9; /* 20 rows. */
  double *x_10; /* 50 rows. */
  double *u_10; /* 20 rows. */
  double *x_11; /* 50 rows. */
  double *u_11; /* 20 rows. */
  double *x_12; /* 50 rows. */
  double *u_12; /* 20 rows. */
  double *x_13; /* 50 rows. */
  double *u_13; /* 20 rows. */
  double *x_14; /* 50 rows. */
  double *u_14; /* 20 rows. */
  double *x_15; /* 50 rows. */
  double *u_15; /* 20 rows. */
  double *x_16; /* 50 rows. */
  double *u_16; /* 20 rows. */
  double *x_17; /* 50 rows. */
  double *u_17; /* 20 rows. */
  double *x_18; /* 50 rows. */
  double *u_18; /* 20 rows. */
  double *x_19; /* 50 rows. */
  double *u_19; /* 20 rows. */
  double *x_20; /* 50 rows. */
  double *u_20; /* 20 rows. */
  double *x_21; /* 50 rows. */
  double *u_21; /* 20 rows. */
  double *x_22; /* 50 rows. */
  double *u_22; /* 20 rows. */
  double *x_23; /* 50 rows. */
  double *u_23; /* 20 rows. */
  double *x_24; /* 50 rows. */
  double *u_24; /* 20 rows. */
  double *x_25; /* 50 rows. */
  double *u_25; /* 20 rows. */
  double *x_26; /* 50 rows. */
  double *u_26; /* 20 rows. */
  double *x_27; /* 50 rows. */
  double *u_27; /* 20 rows. */
  double *x_28; /* 50 rows. */
  double *u_28; /* 20 rows. */
  double *x_29; /* 50 rows. */
  double *u_29; /* 20 rows. */
  double *x_30; /* 50 rows. */
  double *u_30; /* 20 rows. */

  double *u[31];
  double *x[31];
} Vars;

typedef struct Workspace_t {
  double h[1240];
  double s_inv[1240];
  double s_inv_z[1240];
  double b[1500];
  double q[2120];
  double rhs[6100];
  double x[6100];
  double *s;
  double *z;
  double *y;
  double lhs_aff[6100];
  double lhs_cc[6100];
  double buffer[6100];
  double buffer2[6100];

  double K[155220];
  int KKTi[155220];
  int KKTp[6101];
  int P[6100];

  double * D, * Lx;
  int * Lr, *Lp;
  int LNZ;
  double v[6100];

  double gap;
  double optval;

  double ineq_resid_squared;
  double eq_resid_squared;

  double block_33[1];

  /* Pre-op symbols. */
  double quad_645199597568[1];

  int converged;
} Workspace;

typedef struct Settings_t {
  double resid_tol;
  double eps;
  int max_iters;
  int refine_steps;

  int better_start;
  /* Better start obviates the need for s_init and z_init. */
  double s_init;
  double z_init;

  int verbose;
  /* Show extra details of the iterative refinement steps. */
  int verbose_refinement;
  int debug;

  /* For regularization. Minimum value of abs(D_ii) in the kkt D factor. */
  double e;
} Settings;

extern Vars vars;
extern Params p;
extern Workspace w;
extern Settings st;

/* Function definitions in /Users/bodono/Dropbox/research/cvxgen_edit/cvxgen_ldl/qp_solver/ldl.c: */
void fill_KKT(void);
void matrix_multiply(double *r, double *s);
double check_residual(double *target, double *multiplicand);
void ldl_solve(double *target, double *var);

/* Function definitions in /Users/bodono/Dropbox/research/cvxgen_edit/cvxgen_ldl/qp_solver/solver.c: */
double eval_gap(void);
void set_defaults(void);
void setup_pointers(void);
void setup_indexing(void);
void set_start(void);
void fillrhs_aff(void);
void fillrhs_cc(void);
void refine(double *target, double *var);
double calc_ineq_resid_squared(void);
double calc_eq_resid_squared(void);
void better_start(void);
void fillrhs_start(void);
long solve(void);

/* Function definitions in /Users/bodono/Dropbox/research/cvxgen_edit/cvxgen_ldl/qp_solver/matrix_support.c: */
void multbymA(double *lhs, double *rhs);
void multbymAT(double *lhs, double *rhs);
void multbymG(double *lhs, double *rhs);
void multbymGT(double *lhs, double *rhs);
void multbyP(double *lhs, double *rhs);
void fillq(void);
void fillh(void);
void fillb(void);
void pre_ops(void);

/* Function definitions in /Users/bodono/Dropbox/research/cvxgen_edit/cvxgen_ldl/qp_solver/util.c: */
void tic(void);
float toc(void);
float tocq(void);
void printmatrix(char *name, double *A, int m, int n, int sparse);
double unif(double lower, double upper);
float ran1(long*idum, int reset);
float randn_internal(long *idum, int reset);
double randn(void);
void reset_rand(void);

/* Function definitions in /Users/bodono/Dropbox/research/cvxgen_edit/cvxgen_ldl/qp_solver/testsolver.c: */
int main(int argc, char **argv);
void load_default_data(void);

#endif
