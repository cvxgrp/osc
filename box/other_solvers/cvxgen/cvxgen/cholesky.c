#include "cholesky.h"
#include <stdio.h>
#include "string.h"
#include <stdlib.h>

/*
void cholesky_init(int n, int Ap[], int Ar[], int P[], double **info) {
	*info  = (double *) malloc(AMD_INFO * sizeof(double));
	amd_order(n, Ap, Ar, P, (double *) NULL, *info);
}
*/
void cholesky_factor(int n, int Ap[], int Ar[], double Ax[], int **Lp, int **Lr, double **Lx, double **D, int *LNZ) 
{
	*Lp = (int *) malloc((n+1) * sizeof(int));
	int Parent[n], Lnz[n], Flag[n], Pattern[n];
	double Y[n];

	ldl_symbolic(n, Ap, Ar, *Lp, Parent, Lnz, Flag, NULL, NULL);

	*LNZ = *(*Lp + n);

	*D  = (double *) malloc(n    * sizeof(double));
	*Lx = (double *) malloc(*LNZ * sizeof(double));
	*Lr =    (int *) malloc(*LNZ * sizeof(int));

	ldl_numeric(n, Ap, Ar, Ax, *Lp, Parent, Lnz, *Lr, *Lx, *D, Y, Pattern, Flag,NULL,NULL);
}

void cholesky_solve(int n, double *x, double b[], int Lp[], int Lr[], double Lx[], double D[], int P[])
{
	double bp[n];
  int i;
	if (P == NULL) {
		for (i = 0; i < n; i++) { x[i] = b[i]; }
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
