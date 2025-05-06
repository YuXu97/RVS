function [Pp, Qp] = PolyFullGenerators(M,U,V,psi,pc, CM, Cf, indcum)
%U,V: generators of Psi*K*Psi'
%M: model order
%pc: polynomial coefficients
%psi: Psi*zeta
%This function returns the generators of the output kernel of the Wiener
%system, i.e., O = c1^2*(Psi*K*Psi') + c2^2*(Psi*K*Psi').^2 ... + cM^2*(Psi*K*Psi').^M
% plus all the off-diagonals
%pc(i) = a_ici 

[N,r] = size(U);
rk = zeros(M,1);
for m = 1:M
    rk(m) = factorial(r+m-1)/factorial(r-1)/factorial(m);
end

RK = cumsum(rk);

Ua = zeros(r*N,M+1);
Va = zeros(r*N,M+1);

tmpu = ones(N*r,1);
tmpv = ones(N*r,1);
uvec = U(:);
vvec = V(:);
for i = 0:M
    Ua(:,i+1) = tmpu;
    Va(:,i+1) = tmpv;
    if i ~= M
        tmpu = tmpu.*uvec;
        tmpv = tmpv.*vvec;
    end
end

Pp1 = zeros(N,RK(M-1));
Qp1 = zeros(N,RK(M-1));
Pp2 = zeros(N,RK(M-1));
Qp2 = zeros(N,RK(M-1));

PSI = zeros(N,M);
ptmp = ones(N,1);
for i = 1:M
    PSI(:,i) = ptmp;%psi.^(i-1);
    if i~=M
        ptmp = ptmp.*psi;
    end
end

for m = 1:M-1
    [P,Q] = PoweredGenerators(N, r, m, CM, Cf, indcum, Ua, Va);
    eta1 = PSI(:,1:M-m+1)*pc(m:M);
    eta2 = PSI(:,1:M-m)*pc(m+1:M);
    Pp1(:,RK(m)-rk(m)+1:RK(m)) = P.*(eta1*ones(1,rk(m)));
    Qp1(:,RK(m)-rk(m)+1:RK(m)) = Q.*(eta1*ones(1,rk(m)));
    Pp2(:,RK(m)-rk(m)+1:RK(m)) = P.*((eta2.*psi)*ones(1,rk(m)));
    Qp2(:,RK(m)-rk(m)+1:RK(m)) = Q.*((eta2.*psi)*ones(1,rk(m)));
end
[Pend,Qend] = PoweredGenerators(N, r, M, CM, Cf, indcum, Ua, Va);
Pp = [Pp1 -Pp2 pc(M)^2*Pend];
Qp = [Qp1 Qp2 Qend];
end