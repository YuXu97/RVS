clear;
addpath([pwd '/EGRSS']);
addpath([pwd '/BasicFunctions']);
Mn = 2000;
N = 270; 
n = 70;
M = 2;
pc = randn(M,1);
t = (0:Mn-1)';
a = -0.001; b = 0.1; phase = pi/3;
u = exp(a*t).*cos(b*t+phase);
Pi = [exp(a*t).*cos(b*t+phase) exp(a*t).*sin(b*t+phase)];
Rho = [exp(-a*t).*cos(b*t) exp(-a*t).*sin(b*t)];
u = u(end-N+1:end);
Pi = Pi(end-N+1:end,:);
Rho = Rho(1:N,:);

[~,r] = size(Pi);

lam = 0.731857910343507;
U = lam.^(1:N)'; V = ones(N,1);

K0 = tril(U*V')+triu(V*U',1);
K = K0(1:n,1:n);
[~,p] = size(U);



Gb = create_opk([U V],[Pi Rho]);
Ub = Gb(:,1:p+r); Vb = Gb(:,p+r+1:end);

Psi = CalculatePsi(u, n);

O1 = zeros(N-n+1,N-n+1);
for m = 1:M
    O1 = O1 + exp(pc(m))*(Psi*K*Psi').^m;
end
pc = [-0.079022098651632;-6.079230298345516];
[Pp, Qp] = PolyGenerators(M,Ub,Vb,exp(pc));



O2 = tril(Pp*Qp')+triu(Qp*Pp',1);
O2 = O2(n:N,n:N);

% O21 = tril(Pp(n:N,:)*Qp(n:N,:)')+triu(Qp(n:N,:)*Pp(n:N,:)',1);
% norm(O21-O2)


Psi_ir = CalculatePsi_ir(u, n);
O3 = zeros(N-n+1,N-n+1);
for m = 1:M
    O3 = O3 + exp(pc(m))*(Psi_ir*K0*Psi_ir').^m;
end

[Pp1, Qp1] = CalculateOutputKernelGenerators(Pi, Rho, M, 'TC-bd', [pc;lam]);
Pp1 = Pp1(n:N,:);
Qp1 = Qp1(n:N,:);
O = tril(Pp1*Qp1')+triu(Qp1*Pp1',1);
Ot = CalculateOutputKernelValidation(Psi_ir, Psi_ir, M, 'TC-bd', [pc; lam]);

% norm(O1-O2)
% norm(O2-O3)
norm(O2-Ot)

% [Wt,c] = egrss_potrf(Pp',Qp',1);
% L = tril(Pp*Wt,-1) + diag(c);
% 
% logdet_smsp1 = 2*sum(log(abs(diag(L))));

tic
Pp = Pp(n:N,:);
Qp = Qp(n:N,:);
[Q,R] = qr(Pp,0);
bUs = Q;
bVs = Qp*R';
s = sqrt(norm(bUs,'fro')/norm(bVs,'fro'));
bUs = bUs/s;
bVs = bVs*s;
[Wt,c] = egrss_potrf(bUs',bVs',1);
% [Wt1,c1] = egrss_potrf(Pp',Qp',1);
logdet_smsp2 = 2*sum(reallog(c));


y = randn(N-n+1,1);
Liy = egrss_trsv(bUs',Wt,c,y);
nm_Liy_smsp = norm(Liy,2);
toc

tic
Li = chol(O2+eye(N-n+1))';
logdet_direct = 2*sum(log(abs(diag(Li))));
Linv = eye(N-n+1)/Li;
nm_Liy_direct = sqrt(sumsqr(Linv*y));
toc

Oiinvy_smsp = egrss_trsv(bUs',Wt,c,Liy,'T');
Oiinvy_direct = Linv'*Linv*y;
% norm(Oiinvy_smsp-Oiinvy_direct)

l1 = egrss_trsv(bUs',Wt,c,y);
le = egrss_trsv(bUs',Wt,c,ones(N-n+1,1));
h0_smsp = l1'*le/sumsqr(le);

h0_direct = (y'*sum(Linv'*Linv,2))/sum(Linv'*Linv,'all');

% clear;
% N1 = 20; N2 = 20; lam = rand;
% U = lam.^(1:N1)'; V = ones(N2,1);
% K0 = tril(U*V')+triu(V*U',1);
% 
% K1 = zeros(N1,N2);
% for i =1:N1
%     for j =1:N2
%         K1(i,j) = lam^max(i,j);
%     end
% end