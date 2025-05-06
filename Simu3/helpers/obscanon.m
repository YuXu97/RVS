function [sys2, Tinv] = obscanon(sys1)
% OBSCANON  Convert a MISO system to observer canonical form
%   [SYS2, TINV] = OBSCANON(SYS1) converts the MISO system SYS1 to observer
%   canonical form, returned as SYS2. The transformation matrix Tinv is
%   also returned, which is such that z = Tinv*x with z beeing the state
%   vector in SYS2 and x beeing the state vector in SYS1. Furthermore,
%
%     Az = Tinv * Ax * T
%     Cz = Cx   * T
%     Qz = Tinv * Qx * Tinv'
%     Rz = Rx
%
%     with T = inv(Tinv)

if(isa(sys1,'tf'))
    [ny, nu] = size(sys1,1);
    sys1 = ss(sys1); % Convert to SS
else    
    ny = size(sys1.c,1);
    nu = size(sys1.b,2);
end

if(ny ~= 1)
    error('Can only handle MISO systems');
end

% Store the old B-matrix and replace it with identity. This is done to
% obtain the similarity transformation matrix, T.
oldB = sys1.b;
nx = length(sys1.a);
sys1.b = eye(nx);

% Convert to tf
if(isa(sys1,'ss'))
    sys1 = tf(sys1);
end

% Extract polynomials
a = sys1.den{1}'; % The same for all
b = sys1.num;
b = [b{:}];
b = reshape(b,[], nx);

if(any(b(1,:) ~= 0))
    error('Can only handle strictly proper systems');
end

% Make sure that highest order coeff. in denominator is 1
if(a(1) ~= 1)
    a = a/a(1);
    b = b/a(1);
end

% Truncate
a = a(2:end); % Remove initial one
Tinv = b(2:end,:); % Remove initial 0:s - the result is the transformation inverse

% Build observer state space
A = [-a [eye(nx-1) ; zeros(1, nx-1)]];
C = [1 zeros(1, nx-1)];
D = zeros(1,nu);

sys2 = ss(A, Tinv*oldB, C, D, sys1.Ts);
