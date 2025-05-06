clear;

N = 150; n = 100;

t = (0:N-1)';
a = -0.001; b = 0.1; phase = pi/3;
u = exp(a*t).*cos(b*t+phase);
Pi = [exp(a*t).*cos(b*t+phase) exp(a*t).*sin(b*t+phase)];
Rho = [exp(-a*t).*cos(b*t) exp(-a*t).*sin(b*t)];
[~,r] = size(Pi);

lam = 0.5; 
U = lam.^(1:N)'; V = ones(N,1);

K0 = tril(U*V')+triu(V*U',1);
K = K0(1:n,1:n);
[~,p] = size(U);



Gb = create_opk([U V],[Pi Rho]);
Ub = Gb(:,1:p+r); Vb = Gb(:,p+r+1:end);

Psi = CalculatePsi(u, n);

O1 = Psi*K*Psi';
O2 = tril(Ub*Vb')+triu(Vb*Ub',1);
O2 = O2(n:N,n:N);

O3 = Pi(n:N,:)*Rho(1:n,:)'*K*Psi';



norm((O1-O2)./O2)
norm((O1-O3)./O1)
