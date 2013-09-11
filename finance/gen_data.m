%% Generate data for osc_multi_period_portfolio
% data set up (time-invariant):
clear all;
randn('state',0);
rand('state',0);

% set to 1 to run admm in matlab
test_on = 1;

%% Defining the problem's parameters
% n - # of assets
alpha = 1.8;
% n = 10; m = n; T = 30; rho = 0.1;
% n = 30; m = n; T = 60; rho = 0.1;
 n = 50; m = n; T = 100; rho = 0.1;
x_init = zeros(n,1); % initial portfolio
x_term = x_init;

% parameters:
lambda = 0.4; % risk penalty cost
s = 0.2 + 0.1*rand(n,1); % quadratic trading cost
kappa = 0.1*rand(n, 1); % absolute value trading cost

% covariance generation:
sigma_tilde = 0.01*diag(rand(n,1));
% generate correlation matrix, entries vary approx -0.3 to 0.8
temp = randn(n); temp = temp*temp';
l1 = 0; l2 = 10;
v1 = ceil(sprand(n,1,0.8)); v2 = ceil(sprand(n,1,0.8));
Y = temp + l1*(v1*v1') + l2*(v2*v2');
C = diag(1./sqrt(diag(Y)))*Y*diag(1./sqrt(diag(Y)));

%log-covariance
sigma_tilde = C.*(sqrt(diag(sigma_tilde)*ones(1,n)).*sqrt(diag(sigma_tilde)*ones(1,n))');

%log-return
mu = 0.03*rand(n,1);

% mean+covariance
r_bar = exp(mu+0.5*diag(sigma_tilde)); % asset mean returns
sigma = (r_bar*r_bar').*(exp(sigma_tilde)-1); % asset covariance matrix
sigma = (sigma+sigma') / 2; % make sure symmetric, sometimes small error

% matrix sqrt of sigma_tilde (to generate samples)
rt_Sigma_t = sqrtm(sigma_tilde);

%%
clear ee Y C

A = diag(r_bar); B = A;
Q = 2*lambda*sigma;
R = 2*(lambda*sigma+diag(s));
S = 2*lambda*sigma;

mat = [Q S; S' R];

q = zeros(n,1);
r = ones(n,1);
c = zeros(n,1);

RHS = [-repmat([q;r],T+1,1); x_init; repmat(c,T,1)];

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
fprintf(fb,'%6.6f ',kappa);
fprintf(fb,'\n');
fprintf(fb,'%6.6f ',x_term);
fprintf(fb,'\n');
fclose(fb);

%% solve using cvx
tic
cvx_begin
cvx_solver 'sdpt3'
variables X(n,T+1) U(m,T+1)
obj = 0;
for t = 1:T+1
    obj = obj + 0.5 * [X(:,t);U(:,t)]'*mat*[X(:,t);U(:,t)]...
        + q'*X(:,t) + r'*U(:,t) + kappa'*abs(U(:,t));
end
for t=1:T
    X(:,t+1) == A*X(:,t)+B*U(:,t)+c;
end
minimize (obj)
X + U >= 0;
U(:,T+1) + X(:,T+1) == x_term;
X(:,1) == x_init;
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
    [L,D,P] = ldl(M,2e-6);    
    disp('running osc robust state estimation')
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
    tic
    for k = 1:MAX_ITER
        
        %% (x, u) - update
        
        f = reshape(-rho*[x_t+z; u_t+y] + repmat([q;r],1,T+1), (T+1)*(n+m), 1);
        sol = P*(L'\(D\(L\(P'*[-f; h]))));
        
        sol = reshape(sol(1:(T+1)*(n+m)), n+m, T+1);
        
        x = sol(1:n,:);
        u = sol(n+1:end,:);
        
        %% (x_t, u_t) - update
       
        x_t_old = x_t;
        u_t_old = u_t;
		% relaxation
		x_h=x;u_h=u;
        x = alpha*x + (1-alpha)*x_t_old;
        u = alpha*u + (1-alpha)*u_t_old; 
        v = x - z;
        w = u - y; 
        
        x_t = v;
        u_t = pos( 1-(kappa*ones(1,T+1)/rho)./abs(w) ).*w;
        
        u_t2 = pos( 1-(kappa*ones(1,T+1)/(2*rho))./abs((w-v)/2) ).*(w-v)/2;
        x_t2 = -u_t2;
        
        x_t( x_t+u_t<0 ) = x_t2( x_t+u_t<0 );
        u_t( x_t+u_t<0 ) = u_t2( x_t+u_t<0 );
        
        %% (z, y) - update
       
        z = z + x_t - x;
        y = y + u_t - u;
        
        %% diagnostics, reporting, termination checks
		%x=x_h;u=u_h;             
        xu = reshape([x_t; u_t], (T+1)*(n+m),1);  
        objval  = 0.5 * ( xu'*Mat*xu + q_obj'*xu ) + ...
            repmat(kappa,T+1,1)'*abs(reshape(u_t,m*(T+1),1));
   
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
opt = 0;
for t = 1:T+1
    opt = opt + 0.5 * [x_t(:,t);u_t(:,t)]'*mat*[x_t(:,t);u_t(:,t)]...
        + q'*x_t(:,t) + r'*u_t(:,t) + kappa'*abs(u_t(:,t));
end

