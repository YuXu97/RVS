function b = trilv(Ut,Vt,x)
% trilv  Computes matrix-vector product A*x where A is a lower triangular matrix tril(Ut'*Vt).
%
% b = trilv(Ut,Vt,x) computes b = A*x where A is given by tril(Ut'*Vt).
%
% See also: EGRSS_GENERATORS

assert(all(size(Ut) == size(Vt)),'Dimension mismatch: Ut and Vt must be of the same size.')
assert(isvector(x) && length(x) == size(Ut,2),'Dimension mismatch: x must be a vector of length size(Ut,2)')

[p,n] = size(Ut);

b = x;
z = zeros(p,1);
% y = Ut*x;

for k = 1:n
    z = z + Vt(:,k)*b(k);
    %y = y - Ut(:,k)*b(k);
    b(k) = Ut(:,k)'*z;% + Vt(:,k)'*y;
end

end
