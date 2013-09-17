% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize(quad_form(x_0, Q) + quad_form(u_0, R) + quad_form(x_1, Q) + quad_form(u_1, R) + quad_form(x_2, Q) + quad_form(u_2, R) + quad_form(x_3, Q) + quad_form(u_3, R) + quad_form(x_4, Q) + quad_form(u_4, R) + quad_form(x_5, Q) + quad_form(u_5, R) + quad_form(x_6, Q) + quad_form(u_6, R) + quad_form(x_7, Q) + quad_form(u_7, R) + quad_form(x_8, Q) + quad_form(u_8, R) + quad_form(x_9, Q) + quad_form(u_9, R) + quad_form(x_10, Q) + quad_form(u_10, R) + quad_form(x_11, Q) + quad_form(u_11, R) + quad_form(x_12, Q) + quad_form(u_12, R) + quad_form(x_13, Q) + quad_form(u_13, R) + quad_form(x_14, Q) + quad_form(u_14, R) + quad_form(x_15, Q) + quad_form(u_15, R) + quad_form(x_16, Q) + quad_form(u_16, R) + quad_form(x_17, Q) + quad_form(u_17, R) + quad_form(x_18, Q) + quad_form(u_18, R) + quad_form(x_19, Q) + quad_form(u_19, R) + quad_form(x_20, Q) + quad_form(u_20, R) + quad_form(x_21, Q) + quad_form(u_21, R) + quad_form(x_22, Q) + quad_form(u_22, R) + quad_form(x_23, Q) + quad_form(u_23, R) + quad_form(x_24, Q) + quad_form(u_24, R) + quad_form(x_25, Q) + quad_form(u_25, R) + quad_form(x_26, Q) + quad_form(u_26, R) + quad_form(x_27, Q) + quad_form(u_27, R) + quad_form(x_28, Q) + quad_form(u_28, R) + quad_form(x_29, Q) + quad_form(u_29, R) + quad_form(x_30, Q) + quad_form(u_30, R))
%   subject to
%     x_1 == A*x_0 + B*u_0
%     x_2 == A*x_1 + B*u_1
%     x_3 == A*x_2 + B*u_2
%     x_4 == A*x_3 + B*u_3
%     x_5 == A*x_4 + B*u_4
%     x_6 == A*x_5 + B*u_5
%     x_7 == A*x_6 + B*u_6
%     x_8 == A*x_7 + B*u_7
%     x_9 == A*x_8 + B*u_8
%     x_10 == A*x_9 + B*u_9
%     x_11 == A*x_10 + B*u_10
%     x_12 == A*x_11 + B*u_11
%     x_13 == A*x_12 + B*u_12
%     x_14 == A*x_13 + B*u_13
%     x_15 == A*x_14 + B*u_14
%     x_16 == A*x_15 + B*u_15
%     x_17 == A*x_16 + B*u_16
%     x_18 == A*x_17 + B*u_17
%     x_19 == A*x_18 + B*u_18
%     x_20 == A*x_19 + B*u_19
%     x_21 == A*x_20 + B*u_20
%     x_22 == A*x_21 + B*u_21
%     x_23 == A*x_22 + B*u_22
%     x_24 == A*x_23 + B*u_23
%     x_25 == A*x_24 + B*u_24
%     x_26 == A*x_25 + B*u_25
%     x_27 == A*x_26 + B*u_26
%     x_28 == A*x_27 + B*u_27
%     x_29 == A*x_28 + B*u_28
%     x_30 == A*x_29 + B*u_29
%     u_0 <= u_max
%     u_1 <= u_max
%     u_2 <= u_max
%     u_3 <= u_max
%     u_4 <= u_max
%     u_5 <= u_max
%     u_6 <= u_max
%     u_7 <= u_max
%     u_8 <= u_max
%     u_9 <= u_max
%     u_10 <= u_max
%     u_11 <= u_max
%     u_12 <= u_max
%     u_13 <= u_max
%     u_14 <= u_max
%     u_15 <= u_max
%     u_16 <= u_max
%     u_17 <= u_max
%     u_18 <= u_max
%     u_19 <= u_max
%     u_20 <= u_max
%     u_21 <= u_max
%     u_22 <= u_max
%     u_23 <= u_max
%     u_24 <= u_max
%     u_25 <= u_max
%     u_26 <= u_max
%     u_27 <= u_max
%     u_28 <= u_max
%     u_29 <= u_max
%     u_30 <= u_max
%     u_0 >= u_min
%     u_1 >= u_min
%     u_2 >= u_min
%     u_3 >= u_min
%     u_4 >= u_min
%     u_5 >= u_min
%     u_6 >= u_min
%     u_7 >= u_min
%     u_8 >= u_min
%     u_9 >= u_min
%     u_10 >= u_min
%     u_11 >= u_min
%     u_12 >= u_min
%     u_13 >= u_min
%     u_14 >= u_min
%     u_15 >= u_min
%     u_16 >= u_min
%     u_17 >= u_min
%     u_18 >= u_min
%     u_19 >= u_min
%     u_20 >= u_min
%     u_21 >= u_min
%     u_22 >= u_min
%     u_23 >= u_min
%     u_24 >= u_min
%     u_25 >= u_min
%     u_26 >= u_min
%     u_27 >= u_min
%     u_28 >= u_min
%     u_29 >= u_min
%     u_30 >= u_min
%
% with variables
%      u_0  20 x 1
%      u_1  20 x 1
%      u_2  20 x 1
%      u_3  20 x 1
%      u_4  20 x 1
%      u_5  20 x 1
%      u_6  20 x 1
%      u_7  20 x 1
%      u_8  20 x 1
%      u_9  20 x 1
%     u_10  20 x 1
%     u_11  20 x 1
%     u_12  20 x 1
%     u_13  20 x 1
%     u_14  20 x 1
%     u_15  20 x 1
%     u_16  20 x 1
%     u_17  20 x 1
%     u_18  20 x 1
%     u_19  20 x 1
%     u_20  20 x 1
%     u_21  20 x 1
%     u_22  20 x 1
%     u_23  20 x 1
%     u_24  20 x 1
%     u_25  20 x 1
%     u_26  20 x 1
%     u_27  20 x 1
%     u_28  20 x 1
%     u_29  20 x 1
%     u_30  20 x 1
%      x_1  50 x 1
%      x_2  50 x 1
%      x_3  50 x 1
%      x_4  50 x 1
%      x_5  50 x 1
%      x_6  50 x 1
%      x_7  50 x 1
%      x_8  50 x 1
%      x_9  50 x 1
%     x_10  50 x 1
%     x_11  50 x 1
%     x_12  50 x 1
%     x_13  50 x 1
%     x_14  50 x 1
%     x_15  50 x 1
%     x_16  50 x 1
%     x_17  50 x 1
%     x_18  50 x 1
%     x_19  50 x 1
%     x_20  50 x 1
%     x_21  50 x 1
%     x_22  50 x 1
%     x_23  50 x 1
%     x_24  50 x 1
%     x_25  50 x 1
%     x_26  50 x 1
%     x_27  50 x 1
%     x_28  50 x 1
%     x_29  50 x 1
%     x_30  50 x 1
%
% and parameters
%        A  50 x 50
%        B  50 x 20
%        Q  50 x 50   PSD
%        R  20 x 20   PSD
%    u_max   1 x 1    positive
%    u_min   1 x 1    negative
%      x_0  50 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.A, ..., params.x_0, then run
%   [vars, status] = csolve(params, settings)


% Produced by CVXGEN, 2012-08-15 14:12:46 -0700.
% CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2011 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
