// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// header file for cholesky.c

#ifndef CHOLESKY_H_GUARD
#define CHOLESKY_H_GUARD

#include "amd.h"
#include "ldl.h"

void cholesky_init(int n, int Ap[], int Ar[], int P[], double **info);
void cholesky_factor(int n, int Ap[], int Ar[], double Ax[], int P[], int Pinv[],
        int **Lp, int **Lr, double **Lx, double **D, int *LNZ);
void cholesky_solve(int n, double *x, double b[], 
		int Lp[], int Lr[], double Lx[], double D[], int P[]);

#endif
