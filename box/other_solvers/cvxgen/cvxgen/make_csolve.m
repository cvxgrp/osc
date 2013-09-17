% Produced by CVXGEN, 2012-08-15 14:12:46 -0700.
% CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2011 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: make_csolve.m.
% Description: Calls mex to generate the csolve mex file.

%mex -v csolve.c ldl.c matrix_support.c solver.c util.c
mex -c COPTIMFLAGS=-O0 -LLDL/Lib -lldl -ILDL/Include -IUFconfig/ ldl.c
mex COPTIMFLAGS=-Os -LLDL/Lib -lldl -ILDL/Include -IUFconfig/ csolve.c cholesky.c ldl.o matrix_support.c solver.c
