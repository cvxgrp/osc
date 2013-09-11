%% data set up (time-invariant):
clear all;
%randn('seed',0);rand('seed',0);
%n=3;m=2;T=4;
%n=5;m=3;T=5;
n=50;m=20;T=40;

num_steps=1;

noise=0.5*randn(n,num_steps);

A=randn(n);
A=A/max(abs(eig(A)));
B=randn(n,m);
B=1.1*B/max(abs(svd(B)));
c=zeros(n,1);

mat=randn(n+m);mat=mat*mat';
mat(1:n,n+1:end)=zeros(n,m);
mat(n+1:end,1:n)=zeros(m,n);
mat=sparse(mat);

Q=mat(1:n,1:n);
R=mat(n+1:end,n+1:end);
S=mat(1:n,n+1:end);

q=zeros(n,1);
r=zeros(m,1);

x_init=5*randn(n,1);

umax=1;
umin=-1;

alpha=1;

cvx_on=0
x_orig=x_init;

max_its = 5000;

%% fast_mpc 1-step
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
[X,U,telapsed] =...
    fmpc_step(sys,params,zeros(n,T),zeros(m,T),x_orig);
cd ..
%}
%% fast_mpc k-step
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
params.niters=max_its;
params.quiet=true;
params.nsteps=num_steps;

cd ./fast_mpc-0.0.1/
%{
disp('mexing (just in case)')
mex fmpc_step.c
mex fmpc_sim.c
%}
disp('running fast mpc')
tic
[X_hist,U_hist,telapsed] =...
    fmpc_sim(sys,params,[x_init zeros(n,T-1)],zeros(m,T),x_orig,noise);
cd ..
%telapsed
toc

cost_fmpc=0;
for t=1:num_steps
    cost_fmpc=cost_fmpc+X_hist(:,t)'*Q*X_hist(:,t)+U_hist(:,t)'*R*U_hist(:,t);
end
cost_fmpc/num_steps
%}
%% test with cvx:
if cvx_on
    tic
    cvx_begin
    variables X(n,T) U(m,T)
    obj=0;
    for t=1:T
        obj=obj+0.5*[X(:,t);U(:,t)]'*mat*[X(:,t);U(:,t)]...
            +q'*X(:,t)+r'*U(:,t);
    end
    for t=1:T-1
        X(:,t+1)==A*X(:,t)+B*U(:,t)+c;
    end
    minimize(obj)
    U<=umax;
    U>=umin;
    %U(:,T)==0;
    X(:,1)==x_init;
    cvx_end
    toc
end
%% set up matrix:
x_init=x_orig;
rho=10;
%%% factor:
temp1=sparse(mat+rho*eye(n+m));
H=sparse([]);
for t=1:T
    H=blkdiag(H,temp1);
end
temp2=sparse([-A -B]);
temp3=sparse([speye(n) sparse(n,m)]);
M1=sparse([]);M2=sparse([]);
for t=1:T-1
    M1=blkdiag(M1,temp2);
    M2=blkdiag(M2,temp3);
end
M2=blkdiag(M2,temp3);
M1=([[sparse(n,(n+m)*(T-1));M1] sparse(T*n,n+m)]);
M3=M2+M1;
M=[H M3';M3 sparse(n*T,n*T)];
% regularized:
M=[H M3';M3 -1e-6*speye(n*T,n*T)];

%% write matrix to file
cd ../C/
delete dataM;
[i,j,s]=find(M);
i=i-1;
fi = fopen('dataM','w');
fprintf(fi,'%u ',n);fprintf(fi,'%u ',m);fprintf(fi,'%u ',T-1);
fprintf(fi,'%u ',length(i));fprintf(fi,'%6.6f ',rho);fprintf(fi,'%6.6f ',alpha);
fprintf(fi,'\n');
%{
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',j');
fprintf(fi,'\n');
fprintf(fi,'%6.6f ',s');
fprintf(fi,'\n');
fclose(fi);
cd ../matlab/
%}
RHS=[-repmat([q;r],T,1);x_init;repmat(c,T-1,1)];
tmp=full(sum(M~=0));
p=[0 cumsum(tmp)];
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

delete data_prox;
fb = fopen('data_prox','w');
fprintf(fb,'%6.6f ',umax);
fprintf(fb,'%6.6f ',umin);
fclose(fb);
cd ../matlab/

%{
delete dataB;
b=randn((2*n+m)*T,1);
fb = fopen('dataB','w');
fprintf(fb,'%6.6f ',b);
fclose(fb);
cd ../matlab/
%}

%% factorization
%
disp('running factorization')
tic
[L,D,P]=ldl(M,1e-6);
%[L,D,p]=ldl(M,1e-6,'vector');
toc
% tim davis LDL
%{
tic
%p=symamd(M);
%p=1:length(M);
p=cs_amd(M);
[L,D,parent,fl]=ldlsparse(M,p);
L = L + speye (length(L));
toc
%}
inv_D = D\speye(size(D));




%% admm:
disp('running admm mpc')
x_init=x_orig;
tic

w=zeros(m,T);
x=[x_init zeros(n,T-1)];u_t=zeros(m,T);
%x=zeros(n,T);u_t=zeros(m,T);
eps=1e-2;
z=[zeros(T*(n+m),1);x_init;repmat(c,T-1,1)];
x_hist=zeros(n,num_steps);u_hist=zeros(m,num_steps);
for j=1:num_steps
    for k=1:max_its
        x_old=x;
        u_t_old=u_t;
        z(1:T*(n+m))=reshape(rho*[x;w+u_t]-repmat([q;r],1,T),T*(n+m),1);
        sol=P*(L'\(D\(L\(P'*z))));
        
        %sol = cs_ltsolve(L,D\cs_lsolve(L,z(p)));
        %sol = cs_ltsolve(L,inv_D*cs_lsolve(L,z(p)));
        %sol(p)=sol;
        
        %sol=L'\(D\(L\z(p)));
        %sol(p)=sol;
        
        sol=reshape(sol(1:T*(n+m)),n+m,T);
        x=sol(1:n,:);
        u=sol(n+1:end,:);
        
        u_t = max(min(u-w,umax),umin);
        w = w + u_t - u;
        
        rn=norm(u_t-u,'fro');
        dn=rho*(norm([x-x_old;u_t-u_t_old],'fro'));
        if max(rn,dn)<eps
            break
        end
        
    end
    %{
    rn
    dn
    k
    %}
    x_hist(:,j)=x_init;
    u_hist(:,j)=u_t(:,1);
    
    %x_init = A*x_init+B*u_t(:,1)+c+noise(:,j);
    %x = [x(:,2:end) zeros(n,1)];
    %u_t = [u_t(:,2:end) zeros(m,1)];
    %w = [w(:,2:end) zeros(m,1)];
    %z(T*(n+m)+1:T*(n+m)+n) = x_init;
end
toc
%{
cost_admmmpc=0;
for t=1:num_steps
    cost_admmmpc=cost_admmmpc+x_hist(:,t)'*Q*x_hist(:,t)+u_hist(:,t)'*R*u_hist(:,t);
end
cost_admmmpc/num_steps
%}

%% admm (old code):
%{
max_its=1000;
v=zeros(n,T);w=zeros(m,T);
x_t=[x_init zeros(n,T-1)];u_t=zeros(m,T);

eps=1e-2;

%relaxation = 1.6;
tic
for j=1:num_steps
    
    %tic
    x=x_t;u=u_t;
    rnorms=[];dnorms=[];
    
    z=[zeros(T*(n+m),1);x_init;repmat(c,T-1,1)];
    for k=1:max_its
        x_t_old=x_t;
        u_t_old=u_t;
        
        for t=1:T
            z(((t-1)*(n+m)+1):((t-1)*(n+m)+n+m)) = ...
                [rho*(v(:,t)+x_t(:,t))-q;rho*(w(:,t)+u_t(:,t))-r];
        end
        
        sol=P*(L'\(D\(L\(P'*z))));
        sol=reshape(sol(1:T*(n+m)),n+m,T);
        x=sol(1:n,:);
        u=sol(n+1:end,:);
        
        
        x_t = x-v;
        u_t = u-w;
        u_t((u-w)>umax)=umax;
        u_t((u-w)<umin)=umin;
        %u_t(:,T)=zeros(m,1);
        
        
        
        v = v + x_t - x;
        w = w + u_t - u;
        
        
        rn=norm([x_t-x;u_t-u],'fro');
        dn=rho*(norm([x_t-x_t_old;u_t-u_t_old],'fro'));
        %rnorms=[rnorms;rn];
        %dnorms=[dnorms;dn];
        
%{
        if cvx_on
            du=norm(U-u_t,'fro')/norm(U,'fro');
            dx=norm(X-x_t,'fro')/norm(X,'fro');
        end
%}
        if max(rn,dn)<eps
            break
        end
        
    end
    %toc
    %semilogy(rnorms);hold on;semilogy(dnorms,'r');
    %rn
    %dn
    %k
    x_init = A*x_init+B*u_t(:,1)+c+noise(:,j);
    x_t = [x_t(:,2:end) zeros(n,1)];
    u_t = [u_t(:,2:end) zeros(m,1)];
    v = [v(:,2:end) zeros(n,1)];
    w = [w(:,2:end) zeros(m,1)];
    
end
toc
%}
%}