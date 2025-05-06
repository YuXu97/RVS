function [mu, sqS] = mygp(hyp, mean, cov, x, y, xs)
% Exact Gaussian process inference under Gaussian likelihood.

if ischar(mean) || isa(mean, 'function_handle'), mean = {mean}; end  % make cell
if ischar(cov)  || isa(cov,  'function_handle'), cov  = {cov};  end  % make cell

n = size(x,1); % Number of training points

K = feval(cov{:}, hyp.cov, x);         % evaluate covariance matrix
K = (K+K')/2;
m = feval(mean{:}, hyp.mean, x);       % evaluate mean vector (training)
sn2 = exp(2*hyp.lik);                  % noise variance of likGauss

ms = feval(mean{:}, hyp.mean, xs);     % evaluate mean vector (test)
kss = feval(cov{:}, hyp.cov, xs);      % self-covariance'
kss = (kss+kss')/2;
Ks  = feval(cov{:}, hyp.cov, x, xs);   % cross-covariances

% construct a matrix for cholesky
M = [K+sn2*eye(n) Ks ; Ks' kss];
M = M+1e-9*eye(length(M));

R = chol(M);
% [R,err] = chol(M);
% while(err)
%     fprintf('.');
%     lm = min(eig(M));
%     M = M-10*lm*eye(length(M));
%     [R,err] = chol(M);
% end

sqS = R(n+1:end, n+1:end); % square root of covariance is given by the lower right block
A = R(1:n, 1:n);
B = R(1:n, n+1:end);
mu = ms+B'*(A'\(y-m));
% 
% 
% post = infExact(hyp, mean, cov, 'likGauss', x, y);
% 
% 
% 
% alpha = post.alpha; L = post.L; sW = post.sW;
% Ltril = all(all(tril(L,-1)==0));            % is L an upper triangular matrix?
% 
% kss = feval(cov{:}, hyp.cov, xs, xs);  % self-covariance
% Ks  = feval(cov{:}, hyp.cov, x, xs);   % cross-covariances
% ms = feval(mean{:}, hyp.mean, xs);
% mu = ms + Ks'*alpha;                   % predictive means
% 
% if(Ltril)
%     V  = L'\(repmat(sW,1,length(xs)).*Ks);
%     S = kss - V'*V;                       % predictive variances
% else
%     error('L must be upper triangular.');
% end