clear;
n = 80;

lam2 = rand; rho2 = rand; lam3 = rand; rho3 = rand;

[TI,TJ] = meshgrid(1:n, 1:n);
indxi = []; indxj = [];
for i = 0:n-1
    indxi = [indxi; diag(TI,i)];
end
for i = 0:n-1
    indxj = [indxj; diag(TJ,i)];
end
indx = [indxj,indxi];
indxr = indx*[cos(pi/4), -sin(pi/4);sin(pi/4), cos(pi/4)];

tic
P1 = zeros(n*(n+1)/2,n*(n+1)/2);
P2 = zeros(n*(n+1)/2,n*(n+1)/2);
for i = 1:n*(n+1)/2
    for j = i:n*(n+1)/2
        t = indxr(i,:);
        s = indxr(j,:);
        P1(i,j) = lam2^(abs(t(1))+abs(s(1)))*rho2^abs(abs(t(1))-abs(s(1)));
        P2(i,j) = lam3^(abs(t(2))+abs(s(2)))*rho3^abs(abs(t(2))-abs(s(2)));
        P1(j,i) = P1(i,j);
        P2(j,i) = P2(i,j);
    end
end
P = P1.*P2;
toc
% Pk2 = zeros(n*(n+1)/2,n*(n+1)/2);
% for i = 0:n-1
%     fi = (2*n-i+1)*i/2+1;
%     bi = (2*n-i)*(i+1)/2;
%     for j = i:n-1
%         fj = (2*n-j+1)*j/2+1;
%         bj = (2*n-j)*(j+1)/2;
%         Pk2(fi:bi,fj:bj) = lam3^((i+j)*sqrt(2)/2)*rho3^((j-i)*sqrt(2)/2)*ones(bi-fi+1,bj-fj+1);
%     end
% end
% tmp = triu(Pk2,1);
% Pk2 = (tmp+tmp'+diag(diag(Pk2)));

tic
ut2 = (lam3*rho3).^(indxr(:,2));
vt2 = (lam3/rho3).^(indxr(:,2));
Pk2 = tril(ut2*vt2')+triu(vt2*ut2',1);

% ut1 = (lam2*rho2).^(indxr(:,1));
% vt1 = (lam2/rho2).^(indxr(:,1));

v = lam2.^(indxr(:,1));
Pk11 = v*v';

Pk12 = zeros(n*(n+1)/2,n*(n+1)/2);
for i = 1:n*(n+1)/2
    for j = i:n*(n+1)/2
        Pk12(i,j) = rho2^abs(indxr(i,1)-indxr(j,1));
        Pk12(j,i) = Pk12(i,j);
    end
end
Pk1 = Pk11.*Pk12;
Pk = Pk1.*Pk2;
toc
norm(Pk-P)
% t = (1:n)';
% Tr = [t t]*[cos(pi/4), -sin(pi/4);sin(pi/4), cos(pi/4)];
% U1 = (lam2*rho2).^Tr(:,1); U2 = (lam3*rho3).^Tr(:,2);
% V1 = (lam2/rho2).^Tr(:,1); V2 = (lam3/rho3).^Tr(:,2);
% K1 = (tril(U1*V1')+triu(V1*U1',1));
% K2 = (tril(U2*V2')+triu(V2*U2',1));
%
% Pk = zeros(n*(n+1)/2,n*(n+1)/2);
% for i = 0:n-1
%     fi = (2*n-i+1)*i/2+1;
%     bi = (2*n-i)*(i+1)/2;
%     Ui = (lam2*rho2).^indxr(fi:bi,1);
%     Vi = (lam2/rho2).^indxr(fi:bi,1);
%     for j = 0:n-1
%     fj = (2*n-j+1)*j/2+1;
%     bj = (2*n-j)*(j+1)/2;
%     Uj = (lam2*rho2).^indxr(fj:bj,1);
%     Vj = (lam2/rho2).^indxr(fj:bj,1);
%
%     Pk1(fi:bi,fj:bj) = tril(Ui*Vj')+triu(Vi*Uj',1);
%     end
% end


