/* Produced by CVXGEN, 2012-04-22 16:51:14 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */

#include "solver.h"

Vars vars;
Params p;
Workspace w;
Settings st;
#define NUMTESTS 100

void load_data(void);

int main(int argc, char **argv) {
  int num_iters;
#if (NUMTESTS > 0)
  int i;
  double time;
  double time_per;
#endif

  set_defaults();
  setup_indexing();
//  load_default_data();
  load_data();
//  printf("A[1] is %lf\n",params.A[1]);
  /* Solve problem instance for the record. */
  st.verbose = 1;
  num_iters = solve();

#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  st.verbose = 0;

  tic();
  for (i = 0; i < NUMTESTS; i++) {
    solve();
  }
  time = tocq();

  printf("Timed %d solves over %.3f seconds.\n", NUMTESTS, time);
  time_per = time / NUMTESTS;

  if (time_per > 1) {
    printf("Actual time taken per solve: %.3g s.\n", time_per);
  } else if (time_per > 1e-3) {
    printf("Actual time taken per solve: %.3g ms.\n", 1e3*time_per);
  } else {
    printf("Actual time taken per solve: %.3g us.\n", 1e6*time_per);
  }
#endif
#endif

  return 0;
}

void load_data(void) {
	FILE *fp = fopen("data_cvxgen", "r");
	if (fp== NULL){
		printf("Couldn't open file %s\n", "data_cvxgen");
	}
	size_t n,m,T;
	fscanf(fp, "%zu", &n);
	fscanf(fp, "%zu", &m);
	fscanf(fp, "%zu", &T);
	fscanf(fp, "%lf", &p.u_max[0]);
	fscanf(fp, "%lf", &p.u_min[0]);
  int i;
	for(i=0; i<n;i++){
		fscanf(fp, "%lf",&p.x_0[i]);
	}
	//column major order
	for(i=0; i<n*n;i++){
		fscanf(fp, "%lf",&p.Q[i]);
	}
  for(i=0; i<m*m;i++){
		fscanf(fp, "%lf",&p.R[i]);
	}
   for(i=0; i<n*n;i++){
		fscanf(fp, "%lf",&p.A[i]);
	}
	for(i=0; i<n*m;i++){
		fscanf(fp, "%lf",&p.B[i]);
	}
  /*
  printf("Q[0] is %f \n",p.Q[0]);
  printf("R[0] is %f \n",p.R[0]);
	printf("A[0] is %f \n",p.A[0]);
  printf("B[0] is %f \n",p.B[0]);
  printf("u_max is %f \n",p.u_max[0]);
  printf("u_min is %f \n",p.u_min[0]);
  */
}
