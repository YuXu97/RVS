function q = gausspdf(x, Delta, varargin)
% GAUSSPDF  Computes the PDF for a Gaussian
%   Q = GAUSSPDF(X, DELTA) takes a matrix (nx x N) of N column vectors and computes the value of the
%   multivariate (nx-dim) Gaussian PDF with precision DELTA and zero mean.
% 
%   Q = GAUSSPDF(X, DELTA, 'single') takes a matrix (any dim) of  scalars and computes the value of
%   the monovariate Gaussian PDF with precision DELTA and zero mean.
%
%   OBS! The PDF is only computed up to a scalar factor! The same DELTA is used for all vectors.
%
%   Fredrik Lindsten, 2009-10-19
%   lindsten@isy.liu.se


if(nargin == 2)
    % We have N vectors ==> x - nx x N
    % x'*Rinv*x == sum(x.*(Delta*x), 1)
    q = exp(-(1/2)*sum(x.*(Delta*x), 1));
elseif(nargin == 3 && strcmpi(varargin{1}, 'single'))
    % We have N x M x ... scalars
    q = exp(-(1/2)*Delta*x.^2);
end



