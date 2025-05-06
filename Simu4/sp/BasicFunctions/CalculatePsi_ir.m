function [Psi_ir] = CalculatePsi_ir(x, n)
N = length(x);
row = [flip(x(1:n)); zeros(N-n,1)];
col = x(n:N);
Psi_ir = toeplitz(col,row);
end