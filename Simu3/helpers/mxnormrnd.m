function A = mxnormrnd(M,V,K)
[r,c] = size(M);
A = mvnrnd(M(:)', kron(inv(K),V))';
A = reshape(A,[r,c]);