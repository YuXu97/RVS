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