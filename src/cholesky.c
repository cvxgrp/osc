// A Splitting Method for Optimal Control
// by Brendan O'Donoghue, George Stathopoulos and Stephen Boyd
// ===========================================================

// helper functions to perform factorization of KKT matrix
// and subsequent solve steps using the factorization

#include "cholesky.h"
#include <stdio.h>
#include "string.h"
#include <stdlib.h>

void choleskyInit(int n, int Ap[], int Ar[], int P[], double **info) {
	*info  = (double *) malloc(AMD_INFO * sizeof(double));
	amd_order(n, Ap, Ar, P, (double *) NULL, *info);
}

void choleskyFactor(int n, int Ap[], int Ar[], double Ax[], int P[], int Pinv[],
		int **Lp, int **Lr, double **Lx, double **D, int *LNZ) 
{
	*Lp = (int *) malloc((n+1) * sizeof(int));
	int Parent[n], Lnz[n], Flag[n], Pattern[n];
	double Y[n];

	ldl_symbolic(n, Ap, Ar, *Lp, Parent, Lnz, Flag, P, Pinv);

	*LNZ = *(*Lp + n);

	*D  = (double *) malloc(n    * sizeof(double));
	*Lx = (double *) malloc(*LNZ * sizeof(double));
	*Lr =    (int *) malloc(*LNZ * sizeof(int));

	ldl_numeric(n, Ap, Ar, Ax, *Lp, Parent, Lnz, *Lr, *Lx, *D, Y, Pattern, Flag, P, Pinv);
}

void choleskySolve(int n, double *x, double b[], 
		int Lp[], int Lr[], double Lx[], double D[], int P[])
{
	double bp[n];
	if (P == NULL) {
		for (int i = 0; i < n; i++) { x[i] = b[i]; }
		ldl_lsolve(n, x, Lp, Lr, Lx);
		ldl_dsolve(n, x, D);
		ldl_ltsolve(n, x, Lp, Lr, Lx);
	} else {
		ldl_perm(n, bp, b, P);
		ldl_lsolve(n, bp, Lp, Lr, Lx);
		ldl_dsolve(n, bp, D);
		ldl_ltsolve(n, bp, Lp, Lr, Lx);
		ldl_permt(n, x, bp, P);
	}
}
