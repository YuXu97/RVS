function [Pp, Qp] = PolyGenerators(M,U,V,pc, CM, Cf, indcum)
%U,V: generators of Psi*K*Psi'
%M: model order
%pc: polynomial coefficients
%This function returns the generators of the output kernel of the Wiener
%system, i.e., O = a_1^2c1*(Psi*K*Psi') + a_2^2c2*(Psi*K*Psi').^2 ... + a_M^2cM*(Psi*K*Psi').^M
%pc(i) = a_i^2ci >= 0


[N,r] = size(U);
rk = zeros(M,1);
for m = 1:M
    rk(m) = factorial(r+m-1)/factorial(r-1)/factorial(m);
end

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

RK = cumsum(rk);

Pp = zeros(N,RK(M));
Qp = zeros(N,RK(M));
for m = 1:M
    [P,Q] = PoweredGenerators(N, r, m, CM, Cf, indcum, Ua, Va);
    Pp(:,RK(m)-rk(m)+1:RK(m)) = pc(m)*P;
    Qp(:,RK(m)-rk(m)+1:RK(m)) = Q;
end
end