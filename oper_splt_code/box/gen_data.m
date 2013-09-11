%% Generate data for osc_box_control
% data set up (time-invariant):
clear all;
randn('state',0);
rand('state',0);

% set to 1 to run admm in matlab
test_on = 1;

%% Defining the problem's parameters
% System's dimensions and control horizon
alpha = 1.8;
% n - states, m - controls
%n = 5; m = 2; T = 10; rho = 50; x_init = 5*randn(n,1);
n = 20; m = 5; T = 20; rho = 50; x_init = 5*randn(n,1);
%n = 50; m = 20; T = 30; rho = 50; x_init = 5*randn(n,1);

A = randn(n);
A = A/max(abs(eig(A)));
B = randn(n,m);
B = 1.1*B / max(abs(svd(B)));
c = zeros(n,1);

mat = randn(n+m); mat = mat*mat';
mat(1:n,n+1:end) = zeros(n,m);
mat(n+1:end,1:n) = zeros(m,n);
mat = sparse(mat);

Q = mat(1:n,1:n);
R = mat(n+1:end,n+1:end);
S = mat(1:n,n+1:end);

q = zeros(n,1);
r = zeros(m,1);

RHS = [-repmat([q;r],T+1,1); x_init; repmat(c,T,1)];

% prox data
umax = 1;
umin = -1;
%
%% FAST MPC
%{
sys.A=A;
sys.B=B;
sys.Q=Q;
sys.R=R;
sys.xmax=inf*ones(n,1);
sys.xmin=-inf*ones(n,1);
sys.umax=umax*ones(m,1);
sys.umin=umin*ones(m,1);
sys.n=n;
sys.m=m;

params.T=T;
params.Qf=Q;
params.kappa=0.1;
params.niters=5000;
params.quiet=false;

cd ./fast_mpc-0.0.1/
mex fmpc_sim.c -lblas -llapack
mex fmpc_step.c -lblas -llapack

[X,U,telapsed] =...
    fmpc_step(sys,params,[x_init zeros(n,T-1)],zeros(m,T),x_init);
cd ..
disp(sprintf('fast MPC took %f seconds to solve',telapsed))
%}

%% set up matrix - pre-factorize
temp1 = sparse(mat+rho*eye(n+m));
E = sparse([]);
for t = 1:T+1
    E = blkdiag(E,temp1);
end
temp2 = sparse([-A -B]);
temp3 = sparse([speye(n) sparse(n,m)]);
M1 = sparse([]); M2 = sparse([]);
for t = 1:T
    M1 = blkdiag(M1, temp2);
    M2 = blkdiag(M2, temp3);
end
M2 = blkdiag(M2, temp3);
M1 = ( [[sparse(n,(n+m)*T); M1] sparse((T+1)*n,n+m)] );
G = M2 + M1;
% regularized:
M = [ E G'; G -1e-6*speye((T+1)*n,(T+1)*n) ];

[i,j,s] = find(M);
i = i-1;
tmp = full(sum(M~=0));
pw = [0 cumsum(tmp)];

%% write matrix to file
delete data_KKT;
fi = fopen('data_KKT','w');
fprintf(fi,'%u ',n);fprintf(fi,'%u ',m);fprintf(fi,'%u ',T);
fprintf(fi,'%u ',length(i));fprintf(fi,'%6.6f ',rho);fprintf(fi,'%6.6f ',alpha);
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',x_init');
fprintf(fi,'\n');
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',pw);
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',s');
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',RHS');
fprintf(fi,'\n');
fclose(fi);

%% write prox data to file
delete data_prox;
fb = fopen('data_prox','w');
fprintf(fb,'%6.6f ',umax);
fprintf(fb,'%6.6f ',umin);
fclose(fb);

%% solve using cvx
tic
cvx_begin
    cvx_solver 'sdpt3'
    variables X(n,T+1) U(m,T+1)
    obj=0;
    for t=1:T+1
        obj = obj + 0.5 * [X(:,t);U(:,t)]'*mat*[X(:,t);U(:,t)]...
            + q'*X(:,t) + r'*U(:,t);
    end
    for t=1:T
        X(:,t+1) == A*X(:,t) + B*U(:,t) + c;
    end
    minimize (obj)
    subject to 
    X(:,1) == x_init;
    U <= umax;
    U >= umin;
cvx_end
t_cvx = toc;
disp(sprintf('CVX took %f seconds to solve',t_cvx))


%% run admm
if test_on
    Mat = kron(eye(T+1), mat);
    q_obj = kron(ones(T+1,1), 2*[q; r]);
    QUIET    = 0;
    MAX_ITER = 3000;
    EPS_ABS   = 1e-3;
    EPS_REL   = 1e-3;
    disp('running factorization')
    tic
    [L,D,P] = ldl(M,1e-6);    
    disp('running osc box constrained control')
    tic
    % Dimensions
    x = zeros(n,T+1);  % First primal variable - state
    u = zeros(m,T+1);  % First primal variable - input
    x_t = zeros(n,T+1);  % Second primal variable - state
    u_t = zeros(m,T+1);  % Second primal variable - input
    z = zeros(n,T+1);  % Dual variable - state
    y = zeros(m,T+1);  % Dual variable - input
    sol = zeros((n+m)*(T+1),1);  % Vector of KKT solution

    if ~QUIET
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
          'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
    end
    
    % Initializing the KKT right-hand side vector
    f = zeros((n+m)*(T+1),1);
    h = [x_init; repmat(c,T,1)];
    
    for k=1:MAX_ITER
        
        %% (x, u) - update

        f = reshape(-rho*[x; u_t+y],(n+m)*(T+1),1) + ...
                repmat([q; r],T+1,1);

        sol = P*(L'\(D\(L\(P'*[-f; h]))));
        sol = reshape(sol(1:(n+m)*(T+1)),n+m,T+1);

        x = sol(1:n,:);
        u = sol(n+1:end,:);
        
       %% (x_t, u_t) - update
        
       x_t_old = x_t;
       u_t_old = u_t;
       % relaxation
       x = alpha*x + (1-alpha)*x_t_old;
       u = alpha*u + (1-alpha)*u_t_old; 
       v = x - z;
       w = u - y;

       x_t = x;
       u_t = max( min(w,umax), umin );
        
       %% (z, y) - update
       
       z = z + x_t - x;
       y = y + u_t - u;
       
       %% diagnostics, reporting, termination checks    
       xu = reshape([x_t; u_t], (T+1)*(n+m),1);  
       objval  = 0.5 * ( xu'*Mat*xu + q_obj'*xu );
        
       r_norm = norm( [x - x_t; u - u_t], 'fro' );
       s_norm = norm( rho * [x_t - x_t_old; u_t - u_t_old],'fro' );
       eps_pri = sqrt((n+m)*(T+1))*EPS_ABS + ...
                     EPS_REL*max( norm([x_t; u_t ], 'fro'), norm([x; u ], 'fro'));
       eps_dual = sqrt((n+m)*(T+1))*EPS_ABS + EPS_REL*norm([z; y], 'fro');
       if ~QUIET
            fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...q
            r_norm, eps_pri, s_norm, eps_dual, objval);
        end
       
       if (r_norm < eps_pri && s_norm < eps_dual)
             break;
       end   
    end
    toc
end
