function [U,V,Pi] = SI2_kernel(t, p, hyper)
% Simultion Indueced kernel of order p. Computes generator representation
% kernel matrix of simulation induced kernel of order p by a strictly
% monotonic vector t of length n. The expression of the kernel is (41) in
% T. Chen. On kernel design for regularized LTI system identification.
%
% [U,V] = SI2_KERNEL(t,p) returns two matrices U and V of size p-by-n
% such that K = tril(U'*V) + triu(V'*U,1) is the kernel matrix with
% elements
%
%    K(i,j) = sum_{k=0}^{p-1} u_k(t)v_k(s) for t > s

% e^{-\alpha*t(i)}\cos(\beta*t(i))
%
%
% See also: EGRSS_POTRF

% Check inputs
assert(p > 0,'p must be positive')
assert(p == floor(p),'p must be an integer')
assert(isvector(p),'p must be a vector')


% Check monotonicity
if all(diff(t) > 0)
    monotonic = 1;
elseif all(diff(t) < 0)
    monotonic = -1;
else
    error('t must be strictly monotonic')
end

% Convert t to a row vector if it is a column vector
if iscolumn(t)
    t = t';
end

% Preallocated output arrays
U = zeros(p,length(t));
V = zeros(p,length(t));

% Compute coefficients

Pi=zeros(length(t),length(t));

% hyper(1) = exp(-\alpha), log(hyper(1)) = -\alpha, \alpha = \omega_0\xi
% hyper(2) = \xi
% hyper(3) = exp(-\gamma), log(hyper(3)) = -\gamma, b(t) = exp(-\gamma t)
%       hyper(3)
%      --------- = exp(\alpha - \gamma) 
%       hyper(1)
% be = \beta
% theta = \phi

be = -log(hyper(1))/hyper(2)*sqrt(1-hyper(2)^2);
theta = acos(2*log(hyper(3)/hyper(1))/sqrt(4*be^2+4*log(hyper(3)/hyper(1))^2));
for dk = 1:length(t)
    for dj = 1:dk
        Pi(dk,dj) = hyper(1)^(t(dk)+t(dj))*(cos(be*t(dk))-log(hyper(1))/be*sin(be*t(dk)))*(cos(be*t(dj))-log(hyper(1))/be*sin(be*t(dj))) ...
            + hyper(1)^(t(dk)+t(dj))*sin(be*t(dk))*sin(be*t(dj))/be^2 + hyper(1)^(t(dk)+t(dj))*cos(be*(t(dk)-t(dj)))*((hyper(3)/hyper(1))^(2*min(t(dk),t(dj)))-1)/(4*be^2*log(hyper(3)/hyper(1)))...
            + hyper(1)^(t(dk)+t(dj))*(cos(theta+be*(t(dk)+t(dj))) - (hyper(3)/hyper(1))^(2*min(t(dk),t(dj)))*cos(2*be*min(t(dk),t(dj))-theta-be*(t(dk)+t(dj))) )...
            /(2*be^2*sqrt(4*be^2+4*log(hyper(3)/hyper(1))^2));
        Pi(dj,dk) = Pi(dk,dj);
    end
end

% Build U and V
if monotonic == 1
    % t is monotonic increasing
    for k = 0:p-1
        if k==0
        U(k+1,:) = (hyper(1).^t).*cos(be*t);
        V(k+1,:) = (hyper(1).^t).*(cos(be*t+theta) - cos(be*t-theta).*((hyper(3)/hyper(1)).^(2*t)))/(2*be^2*sqrt(4*be^2+4*log(hyper(3)/hyper(1))^2))...                 
                     + (hyper(1).^t).*(-log(hyper(1))/be*sin(be*t) + cos(be*t))...
                     + (hyper(1).^t).*(cos(be*t).*((hyper(3)/hyper(1)).^(2*t)-1)/(4*be^2*log(hyper(3)/hyper(1))));
        end
        if k==1
        U(k+1,:) = (hyper(1).^t).*sin(be*t);
        V(k+1,:) = -(hyper(1).^t).*(sin(be*t+theta) + sin(be*t-theta).*(hyper(3)/hyper(1)).^(2*t))/(2*be^2*sqrt(4*be^2+4*log(hyper(3)/hyper(1))^2))...                 
                     + (hyper(1).^t).*(-log(hyper(1))/be*cos(be*t)+ (log(hyper(1))^2+1)/be^2*sin(be*t))...
                     + (hyper(1).^t).*(sin(be*t).*((hyper(3)/hyper(1)).^(2*t)-1)/(4*be^2*log(hyper(3)/hyper(1))));
        end
                
            
        
    end
else
    % t is monotonic decreasing
    disp('monotonic decreasing case is under construction and only necessary for e.g., the spline kernel');
end

