%% Generate data for osc_supply_chain
% data set up (time-invariant):
clear all;
randn('state',0);
rand('state',0);

% run admm in matlab to verify:
test_on = 1;

% System's dimensions and control horizon
% n - nodes, m - edges
% No. of sources and sinks

alpha = 1.8;
%n = 5; T = 20; numsource = 2; numsink = 2; dist = 0.6; rho=2.5;
%n = 20; T = 20; numsource = 2; numsink = 2; dist = 0.4; rho=2.5;
n = 40; T = 20; numsource = 2; numsink = 2; dist = 0.3; rho=2.5;

%% Constructing a connected graph
% Construct a 2-d random graph placing its nodes inside a unit square.
% The graph has n nodes (warehouses), numsource sources and numsink sinks.
% The m edges are assigned by the algorithm based on a maximum allowable
% distance between warehouses provided by the user in dist.
%
% The constructed graph is tested for connectivity. The user might be asked
% to change dist to acquire a connected graph.


% Checking
while (numsource > n) || (numsink > n)
    error('Check your number of sinks and/or sources; cannot be more than the warehouses');
end

% Build 2-D graph on a unit square
position = rand(2,n);
d = zeros(n);

source_pos = 0.1*rand(2,numsource);
sink_pos = 0.9+0.1*rand(2,numsink);

position=[position source_pos sink_pos];
for i=1:size(position,2)
    for j=1:size(position,2)
        d(i,j) = norm(position(:,i)-position(:,j));
    end
end

% Remove distant (unlike) routes
d(d > dist) = 0;

% Remove links between sources and between sinks
d(n+1:n+numsource, n+1:n+numsource) = 0;
d(n+1+numsource:end, n+1+numsource:end) = 0;

edges = sparse(d);

DG = logical(edges);

% Test for connectivity
for i = 1:size(position,2)
    order{i,:} = graphtraverse(DG,i);
end

for i = 1:size(position,2)
    if (length(order{i,:}) < size(position,2))
        error('Graph is not connected - try to increase maximum allowable distance');
    end
end

%%

% code to plot graph:
%plot(position(1,1:n),position(2,1:n),'.');hold on
%plot(source_pos(1,:),source_pos(2,:),'gx','linewidth',3)
%plot(sink_pos(1,:),sink_pos(2,:),'rx','linewidth',3)

[a,b,c] = find(edges);
for i = 1:length(a)
    %plot([position(1,a(i)),position(1,b(i))],[position(2,a(i)),position(2,b(i))])%,'linewidth',0.01/c(i));
end


% Build the matrices Bminus, Bplus
edges = sparse(d(1:n,1:n));
[row, col, v] = find(edges);
Bminus1 = zeros(n, nnz(edges));
Bplus1 = Bminus1;

for i = 1:nnz(edges)
    Bminus1(row(i),i) = 1;
    Bplus1(col(i),i) = 1;
end


[a_snk, b, d_sink] = find(sparse(d(1:n,n+1+numsource:end)));
sink = zeros(n,length(a_snk));
for i = 1:length(a_snk)
    sink(a_snk(i),i) = 1;
end

[a_src, b, d_source] = find(sparse(d(1:n,n+1:n+numsource)));
source = zeros(n,length(a_src));
for i = 1:length(a_src)
    source(a_src(i),i) = 1;
end
Bplus = [source Bplus1 zeros(n,length(a_snk))];
Bminus = [zeros(n,length(a_src)) Bminus1 sink];
B = Bplus - Bminus;

m = size(B,2);
%disp(sprintf('input dimension m is %i',m));

distances = [d_source; v; d_sink];
% end of graph construction


% Get indices of connected edges to sources and nodes
idx_source = zeros(m,1);
idx_node = zeros(m,1);
for j = 1:m
    if (sum(B(:,j))==1)
        idx_source(j,1) = j;
    end
end
idx_source = nonzeros(idx_source);

% Find the departing edges from each node
for i = 1:n
    idx_depart{i,:} = find(Bminus(i,:)~=0);
end

%% Defining the problem's parameters

% prox data
C = 2;
Ub = 1;

A = eye(n);
c = zeros(n,1);
numinit = round(n/3);

q_t = 0.5*rand(n,1);
q = 0.5*rand(n,1);
% transportation cost scaled according to distances
r = max(distances+0.02*randn(m,1), 0);
% get revenue for edges connected to sinks
r(end-length(a_snk)+1:end) = -15*(0.1*rand(length(a_snk),1)+ones(length(a_snk),1));
% cost of acquisition from sources
r(1:length(a_src)) = 5*(0.1*rand(length(a_src),1)+ones(length(a_src),1));

x_init = C/2*ones(n,1);

RHS = [-repmat([q;r],T+1,1); x_init; repmat(c,T,1)];

%% set up matrix - pre-factorize
temp1=sparse( blkdiag( 2*diag(q_t), zeros(m) ) + rho*eye(n+m) );
E = sparse([]);
for t = 1:T+1
    E = blkdiag(E, temp1);
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
p = [0 cumsum(tmp)];

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
fprintf(fi,'%u ',p);
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',s');
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',RHS');
fprintf(fi,'\n');
fclose(fi);

%% write prox data to file
delete data_prox;
fb = fopen('data_prox','w');
fprintf(fb,'%u ',C); fprintf(fb,'\n');
fprintf(fb,'%u ',Ub); fprintf(fb,'\n');
fprintf(fb,'%u ', length(idx_source));
for i = 1:length(idx_source)
    fprintf(fb,'%u ', idx_source(i)-1);
end
fprintf(fb,'\n');
for i = 1:size(idx_depart,1)
    fprintf(fb,'%u ',length(idx_depart{i}));
    fprintf(fb,'%u ',idx_depart{i}-1);
    fprintf(fb,'\n');
end
fclose(fb);

%% solve using cvx
tic
cvx_begin
cvx_solver 'sdpt3'
variables X(n,T+1) U(m,T+1)
obj = 0;
for t = 1:T+1
    obj = obj + ( q'*X(:,t) + pos(q_t)'*square(X(:,t)) + r'*U(:,t) );
end
minimize ( obj )
subject to
for t = 1:T
    X(:,t+1) == A*X(:,t) + B*U(:,t) + c;
    Bminus*U(:,t) <= X(:,t);
end
X(:,1) == x_init;
Bminus*U(:,T+1) <= X(:,T+1);
X <= C;
U >= 0;
U <= Ub;
cvx_end
t_cvx=toc;
disp(sprintf('CVX took %f seconds to solve',t_cvx))

%% run admm
if test_on
    q_obj = [kron(ones(T+1,1),q); kron(ones(T+1,1),q_t)];
    r_obj = kron(ones(T+1,1),r);
    QUIET    = 0;
    MAX_ITER = 3000;
    EPS_ABS   = 1e-3;
    EPS_REL   = 1e-3;
    disp('running factorization')
    tic
    [L,D,P] = ldl(M,2e-6);
    disp('running osc supply chain')
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
    h = [x_init; zeros(n*T,1)];
    tic
    for k = 1:MAX_ITER
        
        %% (x, u) - update
        
        f = reshape(-rho*[x_t+z; u_t+y],(n+m)*(T+1),1) + ...
            repmat([q; r],T+1,1);
        
        sol = P*(L'\(D\(L\(P'*[-f; h]))));
        sol = reshape(sol(1:(n+m)*(T+1)),n+m,T+1);
        
        x = sol(1:n,:);
        u = sol(n+1:end,:);
        
        %% (x_t, u_t) - update
        
        x_t_old = x_t;
        u_t_old = u_t;
        % relaxation
        %x_h=x;u_h=u;
        x = alpha*x + (1-alpha)*x_t_old;
        u = alpha*u + (1-alpha)*u_t_old;
        v = x - z;
        w = u - y;
        
        % Box constraint on source nodes
        u_t(idx_source,:) = w(idx_source,:);
        for i = 1:n  % For each node..
            wi = w(cell2mat(idx_depart(i))',:);
            vi = v(i,:);
            
            lambda = zeros(1,T+1);
            lambda_old = zeros(1,T+1);
            
            %% Compute lower/upper bound
            ui = max( min(wi, Ub), 0 );
            xi = max( min(vi, C), 0 );
            
            if ( sum((sum(ui,1) > xi)) > 0 )
                lambda(sum(ui,1) > xi) = 1;
                rep_lambda = repmat(lambda,size(wi,1),1);
                ui = max( min(wi - rep_lambda, Ub), 0 );
                xi = max( min(vi + lambda, C), 0 );
                while ( sum((sum(ui,1) > xi)) > 0 )
                    lambda_old(sum(ui,1) > xi) = lambda(sum(ui,1) > xi);
                    lambda(sum(ui,1) > xi) = 2*lambda(sum(ui,1) > xi);
                    rep_lambda = repmat(lambda,size(wi,1),1);
                    ui = max( min(wi - rep_lambda, Ub), 0 );
                    xi = max( min(vi + lambda, C), 0 );
                end
            end
            rep_lambda = repmat(lambda,size(wi,1),1);
            low = lambda_old;
            up = lambda;
            
            
            %% Bisection method
            epsilon = 1e-3;
            
            while (true)
                if ((up-low) <= epsilon)
                    break;
                else
                    lambda((up-low) > epsilon) = (low((up-low) > epsilon) + up((up-low) > epsilon))/2;
                    rep_lambda = repmat(lambda,size(wi,1),1);
                    ui = max( min(wi - rep_lambda, Ub), 0 );
                    xi = max( min(vi + lambda, C), 0 );
                    up(sum(ui,1) <= xi) = lambda(sum(ui,1) <= xi);
                    low(sum(ui,1) > xi) = lambda(sum(ui,1) > xi);
                end
            end
            %  .. compute the state..
            x_t(i,:) = v(i,:) + lambda;
            % ..and the inputs associated with this state.
            u_t(cell2mat(idx_depart(i))',:) = w(cell2mat(idx_depart(i))',:) - rep_lambda;
        end
        
        x_t = min( C, max(x_t, 0) );
        u_t = min( Ub, max(u_t, 0) );
        
        %% (z, y) - update
        
        z = z + x_t - x;
        y = y + u_t - u;
        
        %% diagnostics, reporting, termination checks
        
        %x=x_h;u=u_h;
        
        objval  = q_obj' * [reshape(x_t,n*(T+1),1); square(reshape(x_t,n*(T+1),1))] +...
            r_obj' * reshape(u_t,m*(T+1),1);
        
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
