function y = nonmon(z,d,h)
% Non-monotonic nonlinearity used to test Wiener system identificaiton
% method.

% Coeffs of the polynomial
a = h/(2*d^3);
b = 3*h/(2*d);

% Break points
bp = sqrt((b+1)/(3*a));
c = a*bp^3-b*bp-bp;

ind0 = z < -bp;
ind1 = z >= -bp & z <= bp;
ind2 = z > bp;

y = zeros(1,length(z));
y(ind0) = z(ind0)-c;
y(ind1) = a*z(ind1).^3 - b*z(ind1);
y(ind2) = z(ind2)+c;

y = y';
end