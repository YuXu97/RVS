clear;
tol = sqrt(eps);

lam = 0.9;
rho = 0.6;
ts = 1:80;
U = (lam*rho).^(ts)';
V = (lam/rho).^(ts)';
K = tril(U*V')+triu(V*U',1);
[Uk,Sk,~] = svd(K);

[V,D] = eig(K);
d1 = diag(Sk);
d2 = diag(D);
dsq1 = sqrt(d1);
dsq2 = sqrt(d2);

t1 = Uk*d1/norm(dsq1);
t2 = V*d2/norm(dsq2);
V = V*diag(dsq2);

d3 = dsq2;
d3 = d3/norm(dsq2);
t3 = V*d3; 
% real(sqrt(d));