function [imp] = KautzBasis(tau,xi,n)
% This funtion generates the impusle response (first n elements) of the
% tau-th order kautz basis with the generating pole xi
% Input: tau: the order of the Kautz basis, integer
% xi: the generating pole, xi in (-1,1)
% n: the memory length of the impulse response, integer

% p = xi*ones(tau,1);
% z = [0;ones(tau-1,1)/xi];
% k = sqrt(1-xi^2)*xi^(tau-1);
% sys = zpk(z,p,k,1);
% imp = impulse(sys,0:n-1);

c = zeros(tau,1);
for k = 0:tau-1
    c(k+1) = nchoosek(tau-1,k); 
end
b = sqrt(1-xi^2)*flip(c.*(-xi).^(0:tau-1)');
d = zeros(tau+1,1);
for k = 0:tau
    d(k+1) = nchoosek(tau,k); 
end
a = flip(d.*(-xi).^(tau:-1:0)');
imp = impz(b,a,n);
end