%% Test script for testing egrss_potr(), egrss_trmv(), and egrss_trsv(). 
clear,clc,close all

%% Generate p'th order spline kernel matrix

n = 2000;
p = 2;
t = linspace(0.1,1.0,n);

% Uncomment the two lines below to generate a stable spline kernel matrix
% beta = 0.2;
% t = exp(-beta*t);

[U,V,Pi] = SI2_kernel(t,p,[0.1 0.5 0.2]);

% Compute implicit Cholesky factorization of kernel matrix
W = egrss_potrf(U,V);

%% Form explicit kernel matrix and test Cholesky factorization

% Form explicit kernel matrix 
K = tril(U'*V);
K = K + triu(K',1);

fprintf(1,'n = %i\n',n)
fprintf(1,'p = %i\n',p)
fprintf(1,'cond(K) = %.3e\n',cond(K))

% Compute relative error for creating the kernel matrix
relerr = norm(K-Pi,'fro')/norm(Pi,'fro');
fprintf(1,'\ncreating the matrix  \n')
fprintf(1,'   relative error = %.3e\n',relerr)


% Compute explicit Cholesky factorization of K
tic;
for k=1:2000
Lref = chol(K,'lower');
end
toc
% Compute explicit Cholesky factorization from U and W
tic;
for k=1:2000
L = tril(U'*W);
end
toc

%% Test EGRSS functions 

% Compute relative error for Cholesky factorization
relerr = norm(Lref-L,'fro')/norm(L,'fro');
fprintf(1,'\nCholesky factorization\n')
fprintf(1,'   relative error = %.3e\n',relerr)

% Compute relative error for Cholesky product
Vt = egrss_trtr(U,W);
relerr = norm(Vt-V,'fro')/norm(V,'fro');
fprintf(1,'\nCholesky product\n')
fprintf(1,'   relative error = %.3e\n',relerr)

% Generate random vector
x = randn(n,1);

% Compute y = L*x and compare with reference
y = egrss_trmv(U,W,x);
yref = Lref*x;
fprintf(1,'\nMatrix-vector product (L*x)\n')
fprintf(1,'   relative error = %.3e\n',norm(y-yref)/norm(yref))

% Compute y = L'*x and compare with reference
y = egrss_trmv(U,W,x,'T');
yref = Lref'*x;
fprintf(1,'\nMatrix-vector product (L''*x)\n')
fprintf(1,'   relative error = %.3e\n',norm(y-yref)/norm(yref))

% Compute y = L\x
y = egrss_trsv(U,W,x);
yref = Lref\x;
fprintf(1,'\nForward substitution (L\\x)\n')
fprintf(1,'   relative error = %.3e\n',norm(y-yref)/norm(yref))

% Compute y = L'\x
y = egrss_trsv(U,W,x,'T');
yref = Lref'\x;
fprintf(1,'\nBackward substitution (L''\\x)\n')
fprintf(1,'   relative error = %.3e\n',norm(y-yref)/norm(yref))

% Compute log(det(L))
a = egrss_logdet(U,W);
aref = sum(log(diag(Lref)));
fprintf(1,'\nLog. determinant (log(det(L)))\n')
fprintf(1,'   relative error = %.3e\n',abs(a-aref)/abs(aref))
