% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize(quad_form(x_0, Q) + quad_form(u_0, R) + quad_form(x_1, Q) + quad_form(u_1, R) + quad_form(x_2, Q) + quad_form(u_2, R) + quad_form(x_3, Q) + quad_form(u_3, R) + quad_form(x_4, Q) + quad_form(u_4, R) + quad_form(x_5, Q) + quad_form(u_5, R) + quad_form(x_6, Q) + quad_form(u_6, R) + quad_form(x_7, Q) + quad_form(u_7, R) + quad_form(x_8, Q) + quad_form(u_8, R) + quad_form(x_9, Q) + quad_form(u_9, R) + quad_form(x_10, Q) + quad_form(u_10, R))
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
%
% with variables
%      u_0   2 x 1
%      u_1   2 x 1
%      u_2   2 x 1
%      u_3   2 x 1
%      u_4   2 x 1
%      u_5   2 x 1
%      u_6   2 x 1
%      u_7   2 x 1
%      u_8   2 x 1
%      u_9   2 x 1
%     u_10   2 x 1
%      x_1   5 x 1
%      x_2   5 x 1
%      x_3   5 x 1
%      x_4   5 x 1
%      x_5   5 x 1
%      x_6   5 x 1
%      x_7   5 x 1
%      x_8   5 x 1
%      x_9   5 x 1
%     x_10   5 x 1
%
% and parameters
%        A   5 x 5
%        B   5 x 2
%        Q   5 x 5    PSD
%        R   2 x 2    PSD
%    u_max   1 x 1    positive
%    u_min   1 x 1    positive
%      x_0   5 x 1
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


% Produced by CVXGEN, 2012-04-22 16:51:14 -0700.
% CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2011 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
