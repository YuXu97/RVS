clear;
addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);
addpath([pwd '/EGRSS_mex']);
N = 100; n = 30; M = 3; Nd = N - n + 1; Y = 10*randn(N-n+1,1);

t = (0:N-1)';
a = -0.001; b = 0.1; phase = pi/3;
u = exp(a*t).*cos(b*t+phase);
Pi = [exp(a*t).*cos(b*t+phase) exp(a*t).*sin(b*t+phase)];
Rho = [exp(-a*t).*cos(b*t) exp(-a*t).*sin(b*t)];
[~,r] = size(Pi);

pc = randn(M,1);

lam = 0.8;
U = lam.^(1:N)'; V = ones(N,1);

K0 = tril(U*V')+triu(V*U',1);
K = K0(1:n,1:n);
[~,p] = size(U);

Gb = create_opk([U V],[Pi Rho]);
Ub = Gb(:,1:p+r); Vb = Gb(:,p+r+1:end);

O1 = tril(Ub*Vb') + triu(Vb*Ub',1);
O1 = O1(n:N,n:N);
Psi = CalculatePsi_ir(u, n);


zeta = lam.^(1:n)';
% psi = Psi*zeta;
psi_ir = CalculatePsi_ir(u,n)*(lam.^(1:N)');
XI = psi_ir*ones(1,N-n+1);
% XI1 = repmat(psi_ir,1,N-n+1);

O = zeros(N-n+1,N-n+1);
for i = 1:M
    O = O + pc(i)^2*(Psi*K0*Psi').^i;
end

for i = 1:M
    for j = i+1:M
        TMP = ((Psi*K0*Psi').^i).*(XI.^(j-i));
        O = O + pc(i)*pc(j)*(TMP+TMP');
    end
end

[Pp, Qp] = PolyFullGenerators(M,Ub(n:N,:),Vb(n:N,:),psi_ir,pc);

Osmsp = tril(Pp*Qp')+triu(Qp*Pp',1);

norm(O-Osmsp)

% 
% 
% % zeta = lam.^(1:N)';
% % uu = zeta'*Rho*Pi(n:N,:)';
% % 
% % norm(uu'-psi_ir)
% 
% 
% Oi = O + eye(Nd);
% 
% try
%     L = chol(Oi)';
% catch
%     try
%         L = chol(Oi+eps*eye(Nd))';
%     catch
%         try
%             L = chol(Oi+1e-4*eye(Nd))';
%             chol_flag = 0;
%         catch
%             obj = 1/eps;
%             Oiinv = zeros(Nd,Nd);
%             return
%         end
%     end
% end
% 
% Linv = eye(Nd)/L;
% Oiinv =  Linv'*Linv;
% h0 = (Y'*sum(Oiinv,2))/sum(Oiinv,'all');
% Yn = Y - h0;
% qf = sumsqr(Linv*Yn);
% % qf_t = Y'*Linv'*Linv*Y;
% ld = 2*sum(log(abs(diag(L))));
% obj = Nd*(1-log(Nd))+Nd*log(qf)+ld;
% 
% 
% [Q,R] = qr(Pp,0);
% bUs = Q;
% bVs = Qp*R';
% s = sqrt(norm(bUs,'fro')/norm(bVs,'fro'));
% bUs = bUs/s;
% bVs = bVs*s;
% 
% [Wt,c] = egrss_potrf(bUs',bVs',1);
% logdet = 2*sum(reallog(c));
% 
% l1 = egrss_trsv(bUs',Wt,c,Y);
% le = egrss_trsv(bUs',Wt,c,ones(Nd,1));
% h0_sp = l1'*le/sumsqr(le);
% 
% Ysp = Y - h0_sp;
% 
% Liy = egrss_trsv(bUs',Wt,c,Ysp);
% nm_Liy = norm(Liy,2);
% 
% obj_sp = Nd - Nd*log(Nd) + 2*Nd*log(nm_Liy) + logdet;