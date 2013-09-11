% Generate data for osc_robust_state_estimation
% data set up (time-invariant):
clear all;
randn('state',0);
rand('state',0);

% set to 1 to run admm in matlab
test_on = 1;

%% Defining the problem's parameters
% System's dimensions and control horizon
alpha = 1.8;
% n - states, m - controls, p - outputs
% n = 10; m = n; p = 5; T = 30; rho = 0.1; x_init = 10*randn(n,1);
% n = 30; m = n; p = 10; T = 60; rho = 0.1; x_init = 10*randn(n,1);
 n = 50; m = n; p = 20; T = 100; rho = 0.1; x_init = 10*randn(n,1);

A = 2*randn(n);
A = A / (max(abs(eig(A)))+0.1);
B = eye(m);
C = randn(p,n);
C = C / max( max(C) );
Q = C'*C;
R = zeros(m);
S = zeros(n,m);
c = zeros(n,1);
mat = [Q S; S' R];

x = zeros(n,T+1);
ym = zeros(p,T+1);

% Outliers for the process noise:
u = randn(m,T+1);
prob = 0.25;  % probability with which outliers occur
ts_perturbed = randperm(T+1);
ts_perturbed = sort(ts_perturbed(1:round(prob*(T+1))));
u(:,ts_perturbed) = u(:,ts_perturbed) + 10*randn(m,round(prob*(T+1)));
u_true = u; 

% Sample a Gaussian for measurement noise
sig = 1;
v = sqrt(sig)*randn(p,T+1);

% Generating data
x(:,1) = x_init;
for t = 1:T
    ym(:,t) = C*x(:,t) + v(:,t);
    x(:,t+1) = A*x(:,t) + B*u(:,t);
end
ym(:,T+1) = C*x(:,T+1) + v(:,T+1);
x_true = x;

for t = 1:T+1
    CTym(:,t) = [-C'*ym(:,t); zeros(m,1)];
end

RHS = [-reshape(CTym,(n+m)*(T+1),1); x_init; repmat(c,T,1)];

% prox data
Hub_M = 1;  % halfwidth of the penalty function (see help huber_circ)

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
    M1 = blkdiag(M1,temp2);
    M2 = blkdiag(M2,temp3);
end
M2 = blkdiag(M2,temp3);
M1 = ( [[sparse(n,(n+m)*T); M1] sparse((T+1)*n,n+m)] );
G = M2 + M1;
% regularized:
M=[E G';G -1e-6*speye((T+1)*n,(T+1)*n)];

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
fprintf(fb,'%6.6f ',Hub_M);
fclose(fb);

%% solve using cvx
tic
cvx_begin
cvx_solver 'sdpt3'
    variables X(n,T+1) U(n,T+1)
    obj = 0;
    for t = 1:T+1
        obj = obj + 0.5 * (ym(:,t)-C*X(:,t))'*(ym(:,t)-C*X(:,t)) + 0.5 * huber_circ(U(:,t),[],Hub_M);
    end
    minimize (obj)
    subject to 
    for t = 1:T
        X(:,t+1) == A*X(:,t) + B*U(:,t);
    end
    X(:,1) == x_init;
cvx_end
t_cvx = toc;
disp(sprintf('CVX took %f seconds to solve',t_cvx))


%% run admm
if test_on
    Cobj = kron(eye(T+1),C);
    QUIET    = 0;
    MAX_ITER = 3000;
    EPS_ABS   = 1e-3;
    EPS_REL   = 1e-3;
    disp('running factorization')
    tic
    [L,D,P] = ldl(M,2e-6);    
    disp('running osc robust state estimation')
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
    
    for k = 1:MAX_ITER
        
       %% (x, u) - update
        
       f = reshape(-rho*[x_t+z; u_t+y],(n+m)*(T+1),1) - RHS(1:(n+m)*(T+1));
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

       x_t = v;
       % shrinkage for the inputs 
       for t = 1:T+1
           u_t(:,t) = ( 1 - min( (1/(1+rho)), Hub_M/(rho*norm(w(:,t)))) ) * w(:,t);
       end
        
       %% (z, y) - update
       
       z = z + x_t - x;
       y = y + u_t - u;
        
       %% diagnostics, reporting, termination checks
        
       objval  = 0.5 * (reshape(ym,p*(T+1),1) - Cobj*reshape(x,n*(T+1),1))' * ...
                 (reshape(ym,p*(T+1),1) - Cobj*reshape(x,n*(T+1),1)) + ...
                 0.5 * sum(huber_circ(u,[],Hub_M));
   
       r_norm  = norm( [x - x_t; u - u_t], 'fro' );
       s_norm  = norm( rho * [x_t - x_t_old; u_t - u_t_old], 'fro' );

       eps_pri = sqrt((n+m)*(T+1))*EPS_ABS + ...
                            EPS_REL*max( norm([x_t; u_t ], 'fro'), norm([x; u ], 'fro'));
       eps_dual = sqrt((n+m)*(T+1))*EPS_ABS + EPS_REL*norm([z; y ], 'fro');
        
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
