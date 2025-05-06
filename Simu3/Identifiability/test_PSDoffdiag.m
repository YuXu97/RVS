clear;
E = zeros(1000,1);
for repi = 1:1
N = 200; n =15;
w = mvnrnd(zeros(n,1),eye(n))';
M = 3;
lam = rand; rho = rand; c = 10*rand;
u = rand(N,1);
a = 10*rand(M,1)-5;

Psi = CalculatePsi(u, n);


%xi = Psi*sqrt(c)*((lam*rho).^[0:n-1]'.*cos((0:n-1)'));
xi = Psi*sqrt(c)*(lam*rho).^[1:n]';
%xi = Psi*sqrt(c)*sqrt(1-rho^2)/(1-rho)*((lam.^[0:n-1]').*((1-rho.^[1:n]').*w));
Xi = repmat(xi,1,N-n+1);

U = (lam*rho).^(1:n)';
V = (lam/rho).^(1:n)';
K = c*(tril(U*V')+triu(V*U',1));

% t=lam.^(0:n-1);
% triuP = triu(repmat(t,n,1),1);
% K = c*(triuP + triuP' + diag(t));
Kc = chol(K)';
tmp = Psi*Kc;


O = zeros(N-n+1,N-n+1);


for i = 1:M
    %     O = O + a(i)^2*(xi*xi').^i;
    O = O + a(i)^2*(tmp*tmp').^i;
end


% for i = 1:M
%     for j = i+1:M
%         TMP = (xi.^i)*(xi.^j)';
%         O = O + a(i)*a(j)*(TMP + TMP');
%     end
% end

for i = 1:M
    for j = i+1:M
        TMP = (tmp*tmp').^i.*(Xi.^(j-i));
        O = O + a(i)*a(j)*(TMP + TMP');
    end
end


norm((O-O')./O)
e = eig(O);
min(eig(O))
end




function [Psi] = CalculatePsi(x, n)
%Input: x --> input, n --> memory length
% N= length(x);
% Psi = zeros(N-n,n);
% for j = 1:n
%     Psi(:,j) = x(n-j+1:N-j);
% end
N= length(x);
Psi = zeros(N-n+1,n);

for j=1:n
    Psi(:,j) = x(n-j+1:N-j+1);
end

end