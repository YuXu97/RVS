clear;

N = 200; n = 70; M = 2; Y = randn(N-n+1,1);Nd = N-n+1;
pc = [-1.143281631659069;-0.998304291931689].^3;

% t = (0:N-1)';
% a = -0.001; b = 0.1; phase = pi/3;
% u = exp(a*t).*cos(b*t+phase);
% Pi = [exp(a*t).*cos(b*t+phase) exp(a*t).*sin(b*t+phase)];
% Rho = [exp(-a*t).*cos(b*t) exp(-a*t).*sin(b*t)];
d = load(['Databank/data_N' int2str(10000) '_repi=' int2str(1) '.mat']);
data = d.datainfo.data(1:N,:);
Pi = d.datainfo.Pi(1:N,:);
Rho = d.datainfo.Rho(1:N,:);
u = data(:,1);
[~,r] = size(Pi);

Pi = Pi(n:N,:); Rho = Rho(1:n,:);
Psi = CalculatePsi(u, n); 

lam = 0.893915171593912; rho = 0.999999989620263;
U = (lam*rho).^(1:n)'; V =(lam/rho).^(1:n)';
K = tril(U*V')+triu(V*U',1);

[Ut,S,~] = svd(Rho'*K*Rho);
L = Ut*sqrt(S);
Ub = Pi*L;




% [Pp, ~] = PolyGenerators(M,Ub,Ub,pc);
% O1 = Pp*Pp';
% 
% tmp = Psi*K*Psi';
% O2 = pc(1)*tmp + pc(2)*tmp.^2 + pc(3)*tmp.^3 + pc(4)*tmp.^4;

% norm(O2-O1)

zeta = (lam*rho).^(1:n)';
psi = Pi*Rho'*zeta;
XI = psi*ones(1,N-n+1);
XI1 = repmat(psi,1,N-n+1);
[Pp, Qp] = PolyFullGenerators(M,Ub,Ub,psi,pc);
O1 = Pp*Qp';

O2 = zeros(N-n+1,N-n+1);
for i = 1:M
    O2 = O2 + pc(i)^2*(Psi*K*Psi').^i;
end

for i = 1:M
    for j = i+1:M
        TMP = ((Psi*K*Psi').^i).*(XI.^(j-i));
        O2 = O2 + pc(i)*pc(j)*(TMP+TMP');
    end
end
norm(O1 - O2)


[~,gamma] = size(Pp);
Ig = eye(gamma);
M = Ig+Qp'*Pp;
Mi = Ig/M;

yp = Pp'*Y;
pe = Pp'*ones(Nd,1); qe = Qp'*ones(Nd,1);
h0_sp = (sum(Y)-yp'*Mi*qe)/(Nd-pe'*Mi*qe);
Ysp = Y - h0_sp;
ysp = Pp'*Ysp; ysq = Qp'*Ysp;


qf1 = sumsqr(Ysp) - ysp'*Mi*ysq;

% [Us,Ss,Vs] = svd(eye(gamma)+Qp'*Pp);
% Sinv = Vs*diag(1./diag(Ss))*Us';
% qf2 = sumsqr(Y) - yp'*Sinv*yq;
% norm(eye(gamma)/(eye(gamma)+Qp'*Pp)-Sinv)

% ld1 = logdet(eye(N-n+1)+Pp*Qp');
ld2 = logdet(M);

obj_sp = Nd*(1-log(Nd)) + Nd*log(qf1) + ld2;



Oi2 = O2 + eye(Nd);

try
    L = chol(Oi2)';
catch
    try
        L = chol(Oi2+eps*eye(Nd))';
    catch
        try
            L = chol(Oi2+1e-4*eye(Nd))';
            chol_flag = 0;
        catch
            obj = 1/eps;
            Oiinv = zeros(Nd,Nd);
            return
        end
    end
end

Linv = eye(Nd)/L;
Oiinv =  Linv'*Linv;
h0 = (Y'*sum(Oiinv,2))/sum(Oiinv,'all');
Y = Y - h0;
qf = sumsqr(Linv*Y);
% qf_t = Y'*Linv'*Linv*Y;
ld = 2*sum(log(abs(diag(L))));
obj = Nd*(1-log(Nd))+Nd*log(qf)+ld;





