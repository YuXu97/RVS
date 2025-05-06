function [lb,ub] = lubounds(kernel, M)
 tol = sqrt(eps);


if strcmp(kernel, 'AMLS2os-bd')
    % hp = [c1 c2 ... cM al]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

    lb = [-inf*ones(M,1);  tol;  -tol+1; tol];
    ub = [ inf*ones(M,1);  1-tol; 1-tol; inf];

end

if strcmp(kernel,'SI2od_dc-bd')
    lb = [-inf*ones(M,1); -inf; tol; tol; tol;  tol; tol];
    ub = [inf*ones(M,1); inf; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol];
end

if strcmp(kernel, 'AMLS2od-bd')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

    lb = [-inf*ones(M,1);  tol;  -tol+1; tol; tol];
    ub = [ inf*ones(M,1);  1-tol; 1-tol; inf; inf];

end

if strcmp(kernel,  'SI2od_dc-bd-odd')
    lb = [-inf*ones(M,1); -inf; tol; tol; tol;  tol; tol];
    ub = [inf*ones(M,1); inf; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol];
end

if strcmp(kernel,  'SI2od_dc-bd-odd-polyfix')
    lb = [ -inf; tol; tol; tol;  tol; tol];
    ub = [inf; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol];
end


if strcmp(kernel, 'hDC-bd-odd-polyfix')
     lb = [ -inf;  tol; tol;  tol; tol; tol; tol];
     ub = [inf; 1-tol; 1-tol; 1-tol; inf; inf; inf];
end

if strcmp(kernel, 'hDC_2dc-bd-odd-polyfix')
    lb = [ -inf; -inf;   tol;   tol;   tol;   tol;   tol; tol; tol;  tol];
    ub = [  inf;  inf; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol; inf; inf;  inf];
end

if strcmp(kernel, 'DC-bd')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

    lb = [-inf*ones(M,1);  tol; tol];
    ub = [ inf*ones(M,1);  1-tol; 1-tol];

end

if strcmp(kernel, 'tDC-bd-odd')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

%     lb = [-inf*ones(M,1);  tol; tol];
%     ub = [ inf*ones(M,1);  1-tol; 1-tol];
    lb = [-inf*ones(M,1); -inf;  tol; tol];
    ub = [ inf*ones(M,1); inf; 1-tol; 1-tol];

end

if strcmp(kernel, 'DC-bd-odd')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

%     lb = [-inf*ones(M,1);  tol; tol];
%     ub = [ inf*ones(M,1);  1-tol; 1-tol];

    lb = [-inf*ones(M,1); -inf;  tol; tol];
    ub = [ inf*ones(M,1); inf; 1-tol; 1-tol];

end

if strcmp(kernel, 'SI2od')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

%     lb = [-inf*ones(M,1);  tol; tol];
%     ub = [ inf*ones(M,1);  1-tol; 1-tol];

    lb = [ -inf;  tol; tol; tol];
    ub = [  inf; 1-tol; 1-tol; 1-tol];

end

if strcmp(kernel, '2DC_SI2od') || strcmp(kernel,'2DC_SI2od-bd-odd-polyfix')
    lb = [ -inf; -inf;   tol;   tol;   tol;   tol;   tol;  tol;   tol];
    ub = [  inf;  inf; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol];

end

if strcmp(kernel, 'SI2od_dc') || strcmp(kernel, 'SI2od_tdc')

    lb = [ -inf; tol; tol;   tol; tol; tol];
    ub = [  inf; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol];

end

if strcmp(kernel, 'DC-bd-fm')

     lb = [-inf*ones(M,1)];
     ub = [ inf*ones(M,1)];

end


if strcmp(kernel, 'tDC-bd-odd2')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

     lb = [-inf*ones(M,1);  tol; tol];
     ub = [ inf*ones(M,1);  1-tol; 1-tol];
end


if strcmp(kernel, 'hDC-bd-odd')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

     lb = [-inf*ones(M,1); -inf;  tol; tol;  tol; tol; tol; tol];
     ub = [inf*ones(M,1); inf; 1-tol; 1-tol; 1-tol; inf; inf; inf];

end


if strcmp(kernel, 'hOS')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

     lb = [-inf;   tol; tol; tol; tol];
     ub = [inf; 1-tol; inf; inf; inf];

end


if strcmp(kernel, 'DC-bd-odd-polyfix')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

     lb = [-inf;  tol; tol;  tol; tol; tol; tol];
     ub = [ inf; 1-tol; 1-tol; 1-tol; inf; inf; inf];

end

if strcmp(kernel, 'DC-bd-odd-r1')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

     lb = [tol];
     ub = [inf];

end

if strcmp(kernel, 'DC-bd-odd2')
    % hp = [c1 c2 ... cM lam rho]
    % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
    % cm = am^2 c^(2m)

     lb = [-inf*ones(M,1);  tol; tol];
     ub = [ inf*ones(M,1);  1-tol; 1-tol];

end


if strcmp(kernel, 'DC-opt')
    % hp = [c1 c2 ... cM lam rho]
    % DC-opt: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
    % cm = am c^m
    lb = [tol; -inf*ones(M-1,1);  tol;  tol];
    ub = [inf; inf*ones(M-1,1);  1-tol; 1-tol];
end


if strcmp(kernel, 'DC-dc')
    % hp = [c1 c2 ... cM lam rho]
    % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
    % cm = am c^m
    lb = [tol; -inf*ones(M-1,1);  tol; tol];
    ub = [inf; inf*ones(M-1,1);  1-tol; 1-tol];
end

if strcmp(kernel, 'DC-dcp')
    % hp = [c1 c2 ... cM lam rho]
    % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
    % cm = am c^m
    lb = [tol; -inf*ones(M-1,1);  tol;  tol];
    ub = [inf; inf*ones(M-1,1);  1-tol; 1-tol];
end

if strcmp(kernel, 'DC-dcr')
    % hp = [c1 c2 ... cM lam rho]
    % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
    % cm = am c^m
    lb = [tol; -inf*ones(M-1,1);  tol;  tol; -1+tol];
    ub = [inf; inf*ones(M-1,1);  1-tol; 1-tol; 1-tol];
end

if strcmp(kernel, 'DC-dcs')
    % hp = [c1 c2 ... cM lam rho]
    % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
    % cm = am c^m
    lb = [tol; -inf*ones(M-1,1);  tol;  tol; tol];
    ub = [inf; inf*ones(M-1,1);  1-tol; 1-tol; 1-tol];
end

if strcmp(kernel, 'DC-ob') || strcmp(kernel, 'DC-eig')
    % hp = [c1 c2 ... cM lam rho alpha]

    lb = [tol; -inf*ones(M-1,1);  tol;  tol];
    ub = [inf; inf*ones(M-1,1);  1-tol; 1-tol];

end

if strcmp(kernel, 'DC-ob-ex')
    % hp = [c1 c2 ... cM lam rho alpha]

    lb = [tol; -inf*ones(M-1,1);  tol;  tol];
    ub = [inf; inf*ones(M-1,1);  1-tol; 1-tol];

end

end
