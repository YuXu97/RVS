function [lb,ub] = lubounds(kernel, M)
tol = 1e-3;%sqrt(eps);

if strcmp(kernel, 'DC_mpoly')
    % hp = [c1 c2 ... cM lam1 lam21 lam22 ... lamM1 ...lamMM rho1 rho21 rho22 ... rhoM1 ...rhoMM]
    % number: M^2 + 2M + 4
    lb = [-inf*ones(M,1); tol*ones(M*(M+1)/2,1); tol*ones(M*(M+1)/2,1)];
    ub = [inf*ones(M,1); ones(M*(M+1)/2,1)-tol; ones(M*(M+1)/2,1)];
end

if strcmp(kernel, 'DC_wh_bd') || strcmp(kernel, 'DC_wh_opt') ...
        || strcmp(kernel, 'DC_wh_dc') || strcmp(kernel, 'DC_wh_ob') || strcmp(kernel, 'DC_wh_oba')
    % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2 a0]

    lb = [-inf*ones(M,1); tol; tol; tol; tol;];
    ub = [ inf*ones(M,1); 1-tol; 1-tol; 1; 1];

end

if strcmp(kernel, 'DC_wh_dcs')
    % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2 a0]

    lb = [-inf*ones(M,1); tol; tol; tol; tol; tol];
    ub = [ inf*ones(M,1); 1-tol; 1-tol; 1; 1; 1-tol];

end
end
