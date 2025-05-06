clear;

n = 7;
lam = rand; rho = rand;

[TI,TJ] = meshgrid(1:n, 1:n);
indxi = []; indxj = [];
for i = 0:n-1
    indxi = [indxi; diag(TI,i)];
end
for i = 0:n-1
    indxj = [indxj; diag(TJ,i)];
end
indx = [indxj,indxi];
indx = indx*[cos(pi/4), -sin(pi/4); sin(pi/4), cos(pi/4)];

id1 = indx(1:n,1);
id2 = indx(n+1:2*n-1,1);

P1 = zeros(length(id1),length(id2));
for i = 1:length(id1)
    for j = 1:length(id2)
        P1(i,j) = rho^abs(id1(i)-id2(j));
    end
end

u1 = rho.^id1; v1 = (1/rho).^id1;
u2 = rho.^id2; v2 = (1/rho).^id2;


