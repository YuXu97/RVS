M = 3;

BD = zeros(M*(M-1)/2,(M+2)*(M-1)/2);
for i = 1:M-1
    r1 = 1+i*(i-1)/2;
    r2 = i+i*(i-1)/2;
    c1 = i*(i+1)/2;
    c2 = i*(i+1)/2+i;
    BD(r1:r2,c1:c2) = toeplitz([1;zeros(r2-r1,1)],[1;-1;zeros(c2-c1-1,1)]);
end

A = blkdiag(BD,[zeros(M*(M-1)/2,1),BD]);
A = [zeros(M*(M-1),M+1),A];