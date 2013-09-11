/* Produced by CVXGEN, 2012-03-09 16:42:34 -0800.  */
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
  //load_default_data();
  load_data();
  printf("A[1] is %lf\n",params.A[1]);
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
  params.x_0[5] = 1.58517235573375;
  params.x_0[6] = -1.49765875814465;
  params.x_0[7] = -1.17102848744725;
  params.x_0[8] = -1.79413118679668;
  params.x_0[9] = -0.236760625397454;
  /* Make this a diagonal PSD matrix, even though it's not diagonal. */
  params.Q[0] = 1.02987621087857;
  params.Q[10] = 0;
  params.Q[20] = 0;
  params.Q[30] = 0;
  params.Q[40] = 0;
  params.Q[50] = 0;
  params.Q[60] = 0;
  params.Q[70] = 0;
  params.Q[80] = 0;
  params.Q[90] = 0;
  params.Q[1] = 0;
  params.Q[11] = 1.45683322439471;
  params.Q[21] = 0;
  params.Q[31] = 0;
  params.Q[41] = 0;
  params.Q[51] = 0;
  params.Q[61] = 0;
  params.Q[71] = 0;
  params.Q[81] = 0;
  params.Q[91] = 0;
  params.Q[2] = 0;
  params.Q[12] = 0;
  params.Q[22] = 1.64914404761476;
  params.Q[32] = 0;
  params.Q[42] = 0;
  params.Q[52] = 0;
  params.Q[62] = 0;
  params.Q[72] = 0;
  params.Q[82] = 0;
  params.Q[92] = 0;
  params.Q[3] = 0;
  params.Q[13] = 0;
  params.Q[23] = 0;
  params.Q[33] = 1.27848728264798;
  params.Q[43] = 0;
  params.Q[53] = 0;
  params.Q[63] = 0;
  params.Q[73] = 0;
  params.Q[83] = 0;
  params.Q[93] = 0;
  params.Q[4] = 0;
  params.Q[14] = 0;
  params.Q[24] = 0;
  params.Q[34] = 0;
  params.Q[44] = 1.67625490198013;
  params.Q[54] = 0;
  params.Q[64] = 0;
  params.Q[74] = 0;
  params.Q[84] = 0;
  params.Q[94] = 0;
  params.Q[5] = 0;
  params.Q[15] = 0;
  params.Q[25] = 0;
  params.Q[35] = 0;
  params.Q[45] = 0;
  params.Q[55] = 1.59086281741635;
  params.Q[65] = 0;
  params.Q[75] = 0;
  params.Q[85] = 0;
  params.Q[95] = 0;
  params.Q[6] = 0;
  params.Q[16] = 0;
  params.Q[26] = 0;
  params.Q[36] = 0;
  params.Q[46] = 0;
  params.Q[56] = 0;
  params.Q[66] = 1.02398188237717;
  params.Q[76] = 0;
  params.Q[86] = 0;
  params.Q[96] = 0;
  params.Q[7] = 0;
  params.Q[17] = 0;
  params.Q[27] = 0;
  params.Q[37] = 0;
  params.Q[47] = 0;
  params.Q[57] = 0;
  params.Q[67] = 0;
  params.Q[77] = 1.55885408799088;
  params.Q[87] = 0;
  params.Q[97] = 0;
  params.Q[8] = 0;
  params.Q[18] = 0;
  params.Q[28] = 0;
  params.Q[38] = 0;
  params.Q[48] = 0;
  params.Q[58] = 0;
  params.Q[68] = 0;
  params.Q[78] = 0;
  params.Q[88] = 1.25925244690747;
  params.Q[98] = 0;
  params.Q[9] = 0;
  params.Q[19] = 0;
  params.Q[29] = 0;
  params.Q[39] = 0;
  params.Q[49] = 0;
  params.Q[59] = 0;
  params.Q[69] = 0;
  params.Q[79] = 0;
  params.Q[89] = 0;
  params.Q[99] = 1.41510119701007;
  /* Make this a diagonal PSD matrix, even though it's not diagonal. */
  params.R[0] = 1.28352508177132;
  params.R[4] = 0;
  params.R[8] = 0;
  params.R[12] = 0;
  params.R[1] = 0;
  params.R[5] = 1.693137918313;
  params.R[9] = 0;
  params.R[13] = 0;
  params.R[2] = 0;
  params.R[6] = 0;
  params.R[10] = 1.44045371767074;
  params.R[14] = 0;
  params.R[3] = 0;
  params.R[7] = 0;
  params.R[11] = 0;
  params.R[15] = 1.15686773847496;
  params.A[0] = 0.178596072127379;
  params.A[1] = 1.12125905804547;
  params.A[2] = -0.774545870495281;
  params.A[3] = -1.11216846427127;
  params.A[4] = -0.448114969777405;
  params.A[5] = 1.74553459944172;
  params.A[6] = 1.90398168989174;
  params.A[7] = 0.689534703651255;
  params.A[8] = 1.61133643415359;
  params.A[9] = 1.38300348517272;
  params.A[10] = -0.488023834684443;
  params.A[11] = -1.6311319645131;
  params.A[12] = 0.613643610094145;
  params.A[13] = 0.231363049553804;
  params.A[14] = -0.553740947749688;
  params.A[15] = -1.09978198064067;
  params.A[16] = -0.373920334495006;
  params.A[17] = -0.124239005203324;
  params.A[18] = -0.923057686995755;
  params.A[19] = -0.83282890309827;
  params.A[20] = -0.169254402708088;
  params.A[21] = 1.44213565178771;
  params.A[22] = 0.345011617871286;
  params.A[23] = -0.866048550271161;
  params.A[24] = -0.888089973505595;
  params.A[25] = -0.181511697912213;
  params.A[26] = -1.17835862158005;
  params.A[27] = -1.19448515582771;
  params.A[28] = 0.0561402392697676;
  params.A[29] = -1.65108252487678;
  params.A[30] = -0.0656578705936539;
  params.A[31] = -0.551295150448667;
  params.A[32] = 0.830746487262684;
  params.A[33] = 0.986984892408018;
  params.A[34] = 0.764371687423057;
  params.A[35] = 0.756721655019656;
  params.A[36] = -0.505599503404287;
  params.A[37] = 0.67253921894107;
  params.A[38] = -0.640605344172728;
  params.A[39] = 0.2911754794755;
  params.A[40] = -0.696771367740502;
  params.A[41] = -0.219419802945872;
  params.A[42] = -1.75388427668024;
  params.A[43] = -1.02929831126265;
  params.A[44] = 1.88641042469427;
  params.A[45] = -1.0776631825797;
  params.A[46] = 0.765910043789321;
  params.A[47] = 0.601907432854958;
  params.A[48] = 0.895756557749928;
  params.A[49] = -0.0996455574622748;
  params.A[50] = 0.386655098407451;
  params.A[51] = -1.73212230426869;
  params.A[52] = -1.70975144871107;
  params.A[53] = -1.20409589481169;
  params.A[54] = -1.39255601196584;
  params.A[55] = -1.59958262167422;
  params.A[56] = -1.48282454156458;
  params.A[57] = 0.213110927230614;
  params.A[58] = -1.24874070030449;
  params.A[59] = 1.80840497212483;
  params.A[60] = 0.726447115229707;
  params.A[61] = 0.164078693439085;
  params.A[62] = 0.828722403231591;
  params.A[63] = -0.944453316189946;
  params.A[64] = 1.70690273701491;
  params.A[65] = 1.35677223119988;
  params.A[66] = 0.905277993712149;
  params.A[67] = -0.0790401756583599;
  params.A[68] = 1.36841274350659;
  params.A[69] = 0.979009293697437;
  params.A[70] = 0.64130362559845;
  params.A[71] = 1.65590106802375;
  params.A[72] = 0.534662255150299;
  params.A[73] = -0.536237660589562;
  params.A[74] = 0.211378292601782;
  params.A[75] = -1.21447769319945;
  params.A[76] = -1.23171081442559;
  params.A[77] = 0.902678495731283;
  params.A[78] = 1.13974681372452;
  params.A[79] = 1.88839345473506;
  params.A[80] = 1.40388566816601;
  params.A[81] = 0.174377306383291;
  params.A[82] = -1.64083652190774;
  params.A[83] = -0.0445070215355488;
  params.A[84] = 1.7117453902485;
  params.A[85] = 1.15047279801391;
  params.A[86] = -0.0596230957836474;
  params.A[87] = -0.178882554076455;
  params.A[88] = -1.12805692636259;
  params.A[89] = -1.29114647679271;
  params.A[90] = -1.70550532312257;
  params.A[91] = 1.56957275034837;
  params.A[92] = 0.560706467596236;
  params.A[93] = -1.42667073011471;
  params.A[94] = -0.343492321135171;
  params.A[95] = -1.80356430240851;
  params.A[96] = -1.16250660191055;
  params.A[97] = 0.922832496516153;
  params.A[98] = 0.604491081766398;
  params.A[99] = -0.0840868104920891;
  params.B[0] = -0.900877978017443;
  params.B[1] = 0.608892500264739;
  params.B[2] = 1.82579804526952;
  params.B[3] = -0.257917775299229;
  params.B[4] = -1.71946997964932;
  params.B[5] = -1.76907404870813;
  params.B[6] = -1.66851592480977;
  params.B[7] = 1.83882874901288;
  params.B[8] = 0.163043344745975;
  params.B[9] = 1.34984973067889;
  params.B[10] = -1.31986582305146;
  params.B[11] = -0.958619709084339;
  params.B[12] = 0.767910047491371;
  params.B[13] = 1.58228131256793;
  params.B[14] = -0.637246062159362;
  params.B[15] = -1.74130720803887;
  params.B[16] = 1.45647867764258;
  params.B[17] = -0.836510216682096;
  params.B[18] = 0.96432962559825;
  params.B[19] = -1.36786538119402;
  params.B[20] = 0.779853740563504;
  params.B[21] = 1.36567847612459;
  params.B[22] = 0.908608314986837;
  params.B[23] = -0.563569900546034;
  params.B[24] = 0.906759005960792;
  params.B[25] = -1.44213150327016;
  params.B[26] = -0.744723539067112;
  params.B[27] = -0.321668973268222;
  params.B[28] = 1.50884815577727;
  params.B[29] = -1.38503916571543;
  params.B[30] = 1.52049916099726;
  params.B[31] = 1.19585727688322;
  params.B[32] = 1.88649718831192;
  params.B[33] = -0.529188066786158;
  params.B[34] = -1.18024092436888;
  params.B[35] = -1.0377187186616;
  params.B[36] = 1.31145120568568;
  params.B[37] = 1.86091259437566;
  params.B[38] = 0.795239993521694;
  params.B[39] = -0.0700118329046804;
  params.u_max[0] = 0.574099529362266;
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
