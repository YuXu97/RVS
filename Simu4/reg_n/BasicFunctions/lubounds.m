function [lb,ub] = lubounds(kernel, M)
tol = sqrt(eps);

if strcmp(kernel, 'DC-bd-sp')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

    lb = [-inf*ones(M,1);  tol];
    ub = [ inf*ones(M,1);  1-tol];

end

if strcmp(kernel, 'DC-dc-sp') || strcmp(kernel, 'DC-ob-sp')
    % hp = [c1 c2 ... cM lam rho]
    % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
    % cm = am c^m
    lb = [tol; -inf*ones(M-1,1);  tol];
    ub = [inf; inf*ones(M-1,1);  1-tol];
end

end
