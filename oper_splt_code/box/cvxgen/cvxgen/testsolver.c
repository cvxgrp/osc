/* Produced by CVXGEN, 2012-04-22 16:51:14 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */

#include "solver.h"

Vars vars;
Params params;
Workspace work;
Settings settings;
#define NUMTESTS 1000

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
  settings.verbose = 1;
  num_iters = solve();

#ifndef ZERO_LIBRARY_MODE
#if (NUMTESTS > 0)
  /* Now solve multiple problem instances for timing purposes. */
  settings.verbose = 0;

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

void load_default_data(void) {
  params.x_0[0] = 0.203191610298302;
  params.x_0[1] = 0.832591290472419;
  params.x_0[2] = -0.836381044348223;
  params.x_0[3] = 0.0433104207906521;
  params.x_0[4] = 1.57178781739062;
  /* Make this a diagonal PSD matrix, even though it's not diagonal. */
  params.Q[0] = 1.89629308893344;
  params.Q[5] = 0;
  params.Q[10] = 0;
  params.Q[15] = 0;
  params.Q[20] = 0;
  params.Q[1] = 0;
  params.Q[6] = 1.12558531046384;
  params.Q[11] = 0;
  params.Q[16] = 0;
  params.Q[21] = 0;
  params.Q[2] = 0;
  params.Q[7] = 0;
  params.Q[12] = 1.20724287813819;
  params.Q[17] = 0;
  params.Q[22] = 0;
  params.Q[3] = 0;
  params.Q[8] = 0;
  params.Q[13] = 0;
  params.Q[18] = 1.05146720330083;
  params.Q[23] = 0;
  params.Q[4] = 0;
  params.Q[9] = 0;
  params.Q[14] = 0;
  params.Q[19] = 0;
  params.Q[24] = 1.44080984365064;
  /* Make this a diagonal PSD matrix, even though it's not diagonal. */
  params.R[0] = 1.02987621087857;
  params.R[2] = 0;
  params.R[1] = 0;
  params.R[3] = 1.45683322439471;
  params.A[0] = 0.596576190459043;
  params.A[1] = -0.886050869408099;
  params.A[2] = 0.705019607920525;
  params.A[3] = 0.363451269665403;
  params.A[4] = -1.90407247049134;
  params.A[5] = 0.235416351963528;
  params.A[6] = -0.962990212370138;
  params.A[7] = -0.339595211959721;
  params.A[8] = -0.865899672914725;
  params.A[9] = 0.772551673251985;
  params.A[10] = -0.238185129317042;
  params.A[11] = -1.37252904610015;
  params.A[12] = 0.178596072127379;
  params.A[13] = 1.12125905804547;
  params.A[14] = -0.774545870495281;
  params.A[15] = -1.11216846427127;
  params.A[16] = -0.448114969777405;
  params.A[17] = 1.74553459944172;
  params.A[18] = 1.90398168989174;
  params.A[19] = 0.689534703651255;
  params.A[20] = 1.61133643415359;
  params.A[21] = 1.38300348517272;
  params.A[22] = -0.488023834684443;
  params.A[23] = -1.6311319645131;
  params.A[24] = 0.613643610094145;
  params.B[0] = 0.231363049553804;
  params.B[1] = -0.553740947749688;
  params.B[2] = -1.09978198064067;
  params.B[3] = -0.373920334495006;
  params.B[4] = -0.124239005203324;
  params.B[5] = -0.923057686995755;
  params.B[6] = -0.83282890309827;
  params.B[7] = -0.169254402708088;
  params.B[8] = 1.44213565178771;
  params.B[9] = 0.345011617871286;
  params.u_max[0] = 0.56697572486442;
  params.u_min[0] = 0.555955013247203;
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
	fscanf(fp, "%lf", &params.u_max[0]);
	fscanf(fp, "%lf", &params.u_min[0]);
	int i;
	for(i=0; i<n;i++){
		fscanf(fp, "%lf",&params.x_0[i]);
	}
	//column major order
	for(i=0; i<n*n;i++){
		fscanf(fp, "%lf",&params.Q[i]);
	}
	for(i=0; i<m*m;i++){
		fscanf(fp, "%lf",&params.R[i]);
	}
	for(i=0; i<m*m;i++){
		fscanf(fp, "%lf",&params.A[i]);
	}
	for(i=0; i<m*m;i++){
		fscanf(fp, "%lf",&params.B[i]);
	}
}
