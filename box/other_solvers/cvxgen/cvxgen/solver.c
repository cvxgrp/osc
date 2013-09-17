/* Produced by CVXGEN, 2012-08-15 14:08:14 -0700.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.c. */
/* Description: Main solver file. */

#include "solver.h"

double eval_gap(void) {
  int i;
  double gap;

  gap = 0;
  for (i = 0; i < 1240; i++)
    gap += w.z[i]*w.s[i];

  return gap;
}

void set_defaults(void) {
  st.resid_tol = 1e-6;
  st.eps = 1e-4;
  st.max_iters = 25;
  st.refine_steps = 1;

  st.s_init = 1;
  st.z_init = 1;
  st.debug = 0;
  st.verbose = 1;
  st.verbose_refinement = 0;

  st.better_start = 1;

  st.e = 1e-7; /* regularization amount (1e-7 works well) */
}

void setup_pointers(void) {
  w.y = w.x + 2120;
  w.s = w.x + 3620;
  w.z = w.x + 4860;

  vars.u_0 = w.x + 0;
  vars.u_1 = w.x + 20;
  vars.u_2 = w.x + 40;
  vars.u_3 = w.x + 60;
  vars.u_4 = w.x + 80;
  vars.u_5 = w.x + 100;
  vars.u_6 = w.x + 120;
  vars.u_7 = w.x + 140;
  vars.u_8 = w.x + 160;
  vars.u_9 = w.x + 180;
  vars.u_10 = w.x + 200;
  vars.u_11 = w.x + 220;
  vars.u_12 = w.x + 240;
  vars.u_13 = w.x + 260;
  vars.u_14 = w.x + 280;
  vars.u_15 = w.x + 300;
  vars.u_16 = w.x + 320;
  vars.u_17 = w.x + 340;
  vars.u_18 = w.x + 360;
  vars.u_19 = w.x + 380;
  vars.u_20 = w.x + 400;
  vars.u_21 = w.x + 420;
  vars.u_22 = w.x + 440;
  vars.u_23 = w.x + 460;
  vars.u_24 = w.x + 480;
  vars.u_25 = w.x + 500;
  vars.u_26 = w.x + 520;
  vars.u_27 = w.x + 540;
  vars.u_28 = w.x + 560;
  vars.u_29 = w.x + 580;
  vars.u_30 = w.x + 600;
  vars.x_1 = w.x + 620;
  vars.x_2 = w.x + 670;
  vars.x_3 = w.x + 720;
  vars.x_4 = w.x + 770;
  vars.x_5 = w.x + 820;
  vars.x_6 = w.x + 870;
  vars.x_7 = w.x + 920;
  vars.x_8 = w.x + 970;
  vars.x_9 = w.x + 1020;
  vars.x_10 = w.x + 1070;
  vars.x_11 = w.x + 1120;
  vars.x_12 = w.x + 1170;
  vars.x_13 = w.x + 1220;
  vars.x_14 = w.x + 1270;
  vars.x_15 = w.x + 1320;
  vars.x_16 = w.x + 1370;
  vars.x_17 = w.x + 1420;
  vars.x_18 = w.x + 1470;
  vars.x_19 = w.x + 1520;
  vars.x_20 = w.x + 1570;
  vars.x_21 = w.x + 1620;
  vars.x_22 = w.x + 1670;
  vars.x_23 = w.x + 1720;
  vars.x_24 = w.x + 1770;
  vars.x_25 = w.x + 1820;
  vars.x_26 = w.x + 1870;
  vars.x_27 = w.x + 1920;
  vars.x_28 = w.x + 1970;
  vars.x_29 = w.x + 2020;
  vars.x_30 = w.x + 2070;
}

void setup_indexed_params(void) {
  /* In CVXGEN, you can say */
  /*   parameters */
  /*     A[i] (5,3), i=1..4 */
  /*   end */
  /* This function sets up A[2] to be a pointer to A_2, which is a length-15 */
  /* vector of doubles. */
  /* If you access parameters that you haven't defined in CVXGEN, the result */
  /* is undefined. */

  p.x[0] = p.x_0;
}

void setup_indexed_optvars(void) {
  /* In CVXGEN, you can say */
  /*   variables */
  /*     x[i] (5), i=2..4 */
  /*   end */
  /* This function sets up x[3] to be a pointer to x_3, which is a length-5 */
  /* vector of doubles. */
  /* If you access variables that you haven't defined in CVXGEN, the result */
  /* is undefined. */

  vars.u[0] = vars.u_0;
  vars.x[1] = vars.x_1;
  vars.u[1] = vars.u_1;
  vars.x[2] = vars.x_2;
  vars.u[2] = vars.u_2;
  vars.x[3] = vars.x_3;
  vars.u[3] = vars.u_3;
  vars.x[4] = vars.x_4;
  vars.u[4] = vars.u_4;
  vars.x[5] = vars.x_5;
  vars.u[5] = vars.u_5;
  vars.x[6] = vars.x_6;
  vars.u[6] = vars.u_6;
  vars.x[7] = vars.x_7;
  vars.u[7] = vars.u_7;
  vars.x[8] = vars.x_8;
  vars.u[8] = vars.u_8;
  vars.x[9] = vars.x_9;
  vars.u[9] = vars.u_9;
  vars.x[10] = vars.x_10;
  vars.u[10] = vars.u_10;
  vars.x[11] = vars.x_11;
  vars.u[11] = vars.u_11;
  vars.x[12] = vars.x_12;
  vars.u[12] = vars.u_12;
  vars.x[13] = vars.x_13;
  vars.u[13] = vars.u_13;
  vars.x[14] = vars.x_14;
  vars.u[14] = vars.u_14;
  vars.x[15] = vars.x_15;
  vars.u[15] = vars.u_15;
  vars.x[16] = vars.x_16;
  vars.u[16] = vars.u_16;
  vars.x[17] = vars.x_17;
  vars.u[17] = vars.u_17;
  vars.x[18] = vars.x_18;
  vars.u[18] = vars.u_18;
  vars.x[19] = vars.x_19;
  vars.u[19] = vars.u_19;
  vars.x[20] = vars.x_20;
  vars.u[20] = vars.u_20;
  vars.x[21] = vars.x_21;
  vars.u[21] = vars.u_21;
  vars.x[22] = vars.x_22;
  vars.u[22] = vars.u_22;
  vars.x[23] = vars.x_23;
  vars.u[23] = vars.u_23;
  vars.x[24] = vars.x_24;
  vars.u[24] = vars.u_24;
  vars.x[25] = vars.x_25;
  vars.u[25] = vars.u_25;
  vars.x[26] = vars.x_26;
  vars.u[26] = vars.u_26;
  vars.x[27] = vars.x_27;
  vars.u[27] = vars.u_27;
  vars.x[28] = vars.x_28;
  vars.u[28] = vars.u_28;
  vars.x[29] = vars.x_29;
  vars.u[29] = vars.u_29;
  vars.x[30] = vars.x_30;
  vars.u[30] = vars.u_30;
}

void setup_indexing(void) {
  setup_pointers();
  setup_indexed_params();
  setup_indexed_optvars();
}

void set_start(void) {
  int i;

  for (i = 0; i < 2120; i++)
    w.x[i] = 0;

  for (i = 0; i < 1500; i++)
    w.y[i] = 0;

  for (i = 0; i < 1240; i++)
    w.s[i] = (w.h[i] > 0) ? w.h[i] : st.s_init;

  for (i = 0; i < 1240; i++)
    w.z[i] = st.z_init;
}

double eval_objv(void) {
  int i;
  double objv;

  /* Borrow space in work.rhs. */
  multbyP(w.rhs, w.x);

  objv = 0;
  for (i = 0; i < 2120; i++)
    objv += w.x[i]*w.rhs[i];
  objv *= 0.5;

  for (i = 0; i < 2120; i++)
    objv += w.q[i]*w.x[i];

  objv += w.quad_645199597568[0];

  return objv;
}

void fillrhs_aff(void) {
  int i;
  double *r1, *r2, *r3, *r4;

  r1 = w.rhs;
  r2 = w.rhs + 2120;
  r3 = w.rhs + 3360;
  r4 = w.rhs + 4600;

  /* r1 = -A^Ty - G^Tz - Px - q. */
  multbymAT(r1, w.y);
  multbymGT(w.buffer, w.z);
  for (i = 0; i < 2120; i++)
    r1[i] += w.buffer[i];
  multbyP(w.buffer, w.x);
  for (i = 0; i < 2120; i++)
    r1[i] -= w.buffer[i] + w.q[i];

  /* r2 = -z. */
  for (i = 0; i < 1240; i++)
    r2[i] = -w.z[i];

  /* r3 = -Gx - s + h. */
  multbymG(r3, w.x);
  for (i = 0; i < 1240; i++)
    r3[i] += -w.s[i] + w.h[i];

  /* r4 = -Ax + b. */
  multbymA(r4, w.x);
  for (i = 0; i < 1500; i++)
    r4[i] += w.b[i];
}

void fillrhs_cc(void) {
  int i;

  double *r2;
  double *ds_aff, *dz_aff;

  double mu;
  double alpha;
  double sigma;
  double smu;

  double minval;

  r2 = w.rhs + 2120;
  ds_aff = w.lhs_aff + 2120;
  dz_aff = w.lhs_aff + 3360;

  mu = 0;
  for (i = 0; i < 1240; i++)
    mu += w.s[i]*w.z[i];
  /* Don't finish calculating mu quite yet. */

  /* Find min(min(ds./s), min(dz./z)). */
  minval = 0;
  for (i = 0; i < 1240; i++)
    if (ds_aff[i] < minval*w.s[i])
      minval = ds_aff[i]/w.s[i];
  for (i = 0; i < 1240; i++)
    if (dz_aff[i] < minval*w.z[i])
      minval = dz_aff[i]/w.z[i];

  /* Find alpha. */
  if (-1 < minval)
      alpha = 1;
  else
      alpha = -1/minval;

  sigma = 0;
  for (i = 0; i < 1240; i++)
    sigma += (w.s[i] + alpha*ds_aff[i])*
      (w.z[i] + alpha*dz_aff[i]);
  sigma /= mu;
  sigma = sigma*sigma*sigma;

  /* Finish calculating mu now. */
  mu *= 0.0008064516129032258;

  smu = sigma*mu;

  /* Fill-in the rhs. */
  for (i = 0; i < 2120; i++)
    w.rhs[i] = 0;
  for (i = 3360; i < 6100; i++)
    w.rhs[i] = 0;

  for (i = 0; i < 1240; i++)
    r2[i] = w.s_inv[i]*(smu - ds_aff[i]*dz_aff[i]);
}

void refine(double *target, double *var) {
  int i, j;

  double *residual = w.buffer;
  double norm2;
  double *new_var = w.buffer2;
  for (j = 0; j < st.refine_steps; j++) {
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 6100; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }

#ifndef ZERO_LIBRARY_MODE
    if (st.verbose_refinement) {
      if (j == 0)
        printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
      else
        printf("After refinement we get squared norm %.6g.\n", norm2);
    }
#endif

    /* Solve to find new_var = KKT \ (target - A*var). */
    ldl_solve(residual, new_var);

    /* Update var += new_var, or var += KKT \ (target - A*var). */
    for (i = 0; i < 6100; i++) {
      var[i] -= new_var[i];
    }
  }

#ifndef ZERO_LIBRARY_MODE
  if (st.verbose_refinement) {
    /* Check the residual once more, but only if we're reporting it, since */
    /* it's expensive. */
    norm2 = 0;
    matrix_multiply(residual, var);
    for (i = 0; i < 6100; i++) {
      residual[i] = residual[i] - target[i];
      norm2 += residual[i]*residual[i];
    }

    if (j == 0)
      printf("Initial residual before refinement has norm squared %.6g.\n", norm2);
    else
      printf("After refinement we get squared norm %.6g.\n", norm2);
  }
#endif
}

double calc_ineq_resid_squared(void) {
  /* Calculates the norm ||-Gx - s + h||. */
  double norm2_squared;
  int i;

  /* Find -Gx. */
  multbymG(w.buffer, w.x);

  /* Add -s + h. */
  for (i = 0; i < 1240; i++)
    w.buffer[i] += -w.s[i] + w.h[i];

  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 1240; i++)
    norm2_squared += w.buffer[i]*w.buffer[i];

  return norm2_squared;
}

double calc_eq_resid_squared(void) {
  /* Calculates the norm ||-Ax + b||. */
  double norm2_squared;
  int i;

  /* Find -Ax. */
  multbymA(w.buffer, w.x);

  /* Add +b. */
  for (i = 0; i < 1500; i++)
    w.buffer[i] += w.b[i];

  /* Now find the squared norm. */
  norm2_squared = 0;
  for (i = 0; i < 1500; i++)
    norm2_squared += w.buffer[i]*w.buffer[i];

  return norm2_squared;
}

void better_start(void) {
  /* Calculates a better starting point, using a similar approach to CVXOPT. */
  /* Not yet speed optimized. */
  int i;
  double *x, *s, *z, *y;
  double alpha;

  w.block_33[0] = -1;
  /* Make sure sinvz is 1 to make hijacked KKT system ok. */
  for (i = 0; i < 1240; i++)
    w.s_inv_z[i] = 1;
  fill_KKT();
	fill_colcom();
  ldl_factor();

  fillrhs_start();
  /* Borrow work.lhs_aff for the solution. */
  ldl_solve(w.rhs, w.lhs_aff);
  /* Don't do any refinement for now. Precision doesn't matter too much. */

  x = w.lhs_aff;
  s = w.lhs_aff + 2120;
  z = w.lhs_aff + 3360;
  y = w.lhs_aff + 4600;

  /* Just set x and y as is. */
  for (i = 0; i < 2120; i++)
    w.x[i] = x[i];

  for (i = 0; i < 1500; i++)
    w.y[i] = y[i];

  /* Now complete the initialization. Start with s. */
  /* Must have alpha > max(z). */
  alpha = -1e99;
  for (i = 0; i < 1240; i++)
    if (alpha < z[i])
      alpha = z[i];

  if (alpha < 0) {
    for (i = 0; i < 1240; i++)
      w.s[i] = -z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 1240; i++)
      w.s[i] = -z[i] + alpha;
  }

  /* Now initialize z. */
  /* Now must have alpha > max(-z). */
  alpha = -1e99;
  for (i = 0; i < 1240; i++)
    if (alpha < -z[i])
      alpha = -z[i];

  if (alpha < 0) {
    for (i = 0; i < 1240; i++)
      w.z[i] = z[i];
  } else {
    alpha += 1;
    for (i = 0; i < 1240; i++)
      w.z[i] = z[i] + alpha;
  }
}

void fillrhs_start(void) {
  /* Fill rhs with (-q, 0, h, b). */
  int i;
  double *r1, *r2, *r3, *r4;

  r1 = w.rhs;
  r2 = w.rhs + 2120;
  r3 = w.rhs + 3360;
  r4 = w.rhs + 4600;

  for (i = 0; i < 2120; i++)
    r1[i] = -w.q[i];

  for (i = 0; i < 1240; i++)
    r2[i] = 0;

  for (i = 0; i < 1240; i++)
    r3[i] = w.h[i];

  for (i = 0; i < 1500; i++)
    r4[i] = w.b[i];
}

long solve(void) {
  int i;
  int iter;

  double *dx, *ds, *dy, *dz;
  double minval;
  double alpha;
  w.converged = 0;

  setup_pointers();
  pre_ops();

#ifndef ZERO_LIBRARY_MODE
  if (st.verbose)
    printf("iter     objv        gap       |Ax-b|    |Gx+s-h|    step\n");
#endif

  fillq();
  fillh();
  fillb();

  if (st.better_start)
    better_start();
  else
    set_start();

  for (iter = 0; iter < st.max_iters; iter++) {
    for (i = 0; i < 1240; i++) {
      w.s_inv[i] = 1.0 / w.s[i];
      w.s_inv_z[i] = w.s_inv[i]*w.z[i];
    }

    w.block_33[0] = 0;
    fill_KKT();
    ldl_factor();
    /* Affine scaling directions. */
    fillrhs_aff();
    ldl_solve(w.rhs, w.lhs_aff);
    refine(w.rhs, w.lhs_aff);
    /* Centering plus corrector directions. */
    fillrhs_cc();
    ldl_solve(w.rhs, w.lhs_cc);
    refine(w.rhs, w.lhs_cc);

    /* Add the two together and store in aff. */
    for (i = 0; i < 6100; i++)
      w.lhs_aff[i] += w.lhs_cc[i];

    /* Rename aff to reflect its new meaning. */
    dx = w.lhs_aff;
    ds = w.lhs_aff + 2120;
    dz = w.lhs_aff + 3360;
    dy = w.lhs_aff + 4600;

    /* Find min(min(ds./s), min(dz./z)). */
    minval = 0;
    for (i = 0; i < 1240; i++)
      if (ds[i] < minval*w.s[i])
        minval = ds[i]/w.s[i];
    for (i = 0; i < 1240; i++)
      if (dz[i] < minval*w.z[i])
        minval = dz[i]/w.z[i];

    /* Find alpha. */
    if (-0.99 < minval)
      alpha = 1;
    else
      alpha = -0.99/minval;

    /* Update the primal and dual variables. */
    for (i = 0; i < 2120; i++)
      w.x[i] += alpha*dx[i];
    for (i = 0; i < 1240; i++)
      w.s[i] += alpha*ds[i];
    for (i = 0; i < 1240; i++)
      w.z[i] += alpha*dz[i];
    for (i = 0; i < 1500; i++)
      w.y[i] += alpha*dy[i];

    w.gap = eval_gap();
    w.eq_resid_squared = calc_eq_resid_squared();
    w.ineq_resid_squared = calc_ineq_resid_squared();
#ifndef ZERO_LIBRARY_MODE
    if (st.verbose) {
      w.optval = eval_objv();
      printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e  % 6.4f\n",
          iter+1, w.optval, w.gap, sqrt(w.eq_resid_squared),
          sqrt(w.ineq_resid_squared), alpha);
    }
#endif

    /* Test termination conditions. Requires optimality, and satisfied */
    /* constraints. */
    if (   (w.gap < st.eps)
        && (w.eq_resid_squared <= st.resid_tol*st.resid_tol)
        && (w.ineq_resid_squared <= st.resid_tol*st.resid_tol)
       ) {

      w.converged = 1;
      w.optval = eval_objv();
      return iter+1;
    }
  }

  return iter;
}
