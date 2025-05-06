function [Psi] = CalculatePsi_delay(x, n)
%Input: x --> input, n --> memory length
% N= length(x);
% Psi = zeros(N-n,n);
% for j = 1:n
%     Psi(:,j) = x(n-j+1:N-j);
% end
N= length(x);
Psi = zeros(N-n+1,n-1);

for j=2:n
    Psi(:,j-1) = x(n-j+1:N-j+1);
end

end