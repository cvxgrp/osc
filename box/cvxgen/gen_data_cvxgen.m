%% Generate data for osc_box_control
% data set up (time-invariant):
clear all;
randn('state',0);
rand('state',0);

% set to 1 to run admm in matlab

%% Defining the problem's parameters
% System's dimensions and control horizon
% n - states, m - controls
n = 5; m = 2; T = 10; x_init = 5*randn(n,1);
%n = 20; m = 5; T = 20; rho = 50; x_init = 5*randn(n,1);
%n = 50; m = 20; T = 30; rho = 50; x_init = 5*randn(n,1);

A = randn(n);
A = A/max(abs(eig(A)));
params.A=A;
B = randn(n,m);
B = 1.1*B / max(abs(svd(B)));
params.B=B;

mat = randn(n+m); mat = mat*mat';
mat(1:n,n+1:end) = zeros(n,m);
mat(n+1:end,1:n) = zeros(m,n);
%mat = sparse(mat);

params.Q = mat(1:n,1:n);
params.R = mat(n+1:end,n+1:end);

% prox data
params.umax = 1;
params.umin = -1;

params.x0=x_init;


%% write matrix to file
cd cvxgen
delete data_cvxgen;
fi = fopen('data_cvxgen','w');
fprintf(fi,'%u ',n);
fprintf(fi,'%u ',m);
fprintf(fi,'%u ',T);
fprintf(fi,'%6.6f ',params.umax);
fprintf(fi,'%6.6f ',params.umin);

fprintf(fi,'\n');
fprintf(fi,'%6.6f ',params.x0);
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',params.Q);
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',params.R);
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',params.A);
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',params.B);
fclose(fi);
cd ..
%%
% Using CVXGEN to solve
%urlwrite('http://cvxgen.stanford.edu/download/cvxgen.m', 'cvxgen.m');
%cvxgen(934909444096)

%tic
%[vars, status] = csolve(params);
%toc
