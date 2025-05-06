 function [lb,ub] = lubounds(kernel, M)
tol = sqrt(eps);
if strcmp(kernel,'SI2od_dc-bd')
    lb = [-inf*ones(M,1); -inf; tol; tol; tol;  tol; tol; -inf; tol];
    ub = [inf*ones(M,1); inf; 1-tol; 1-tol; 1-tol; 1-tol; 1-tol; inf; inf];
end

if strcmp(kernel, 'TC_bd')
%     % hp = [a1 a2 ... aM c lam rho sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol; tol;  tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;   -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol;  inf; inf];

    lb = [-inf*ones(M,1); -inf; tol;  -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; inf];
end

if strcmp(kernel, 'TC_opt+') || strcmp(kernel, 'TC_opt-')
%     % hp = [a1 a2 ... aM c lam rho sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol; tol;  tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;   -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol;  inf; inf];

    lb = [-inf*ones(M,1); -inf; tol;  -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; inf];
end


if strcmp(kernel, 'TC_tc+') || strcmp(kernel, 'TC_tc-')
%     % hp = [a1 a2 ... aM c lam rho  sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol; tol;  tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;   -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol;  inf; inf];

    lb = [-inf*ones(M,1); -inf; tol;  -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; inf];
end

if strcmp(kernel, 'DC_bd')
%     % hp = [a1 a2 ... aM c lam rho sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol; tol;  tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;  -1; -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf; inf];

    lb = [-inf*ones(M,1); -inf; tol;  -1; -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf];


end


if strcmp(kernel, 'DC_opt+') || strcmp(kernel, 'DC_opt-')
%     % hp = [a1 a2 ... aM c lam rho sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol;  tol; tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;  -1; -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf; inf];
    lb = [-inf*ones(M,1); -inf; tol;  -1; -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf];

end


if strcmp(kernel, 'DC_dc+') || strcmp(kernel, 'DC_dc-')
%     % hp = [a1 a2 ... aM c lam rho sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol;  tol; tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;  -1; -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf; inf];
    lb = [-inf*ones(M,1); -inf; tol;  -1; -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf];

end


if strcmp(kernel, 'DC_eig+') || strcmp(kernel, 'DC_eig-')
%     % hp = [a1 a2 ... aM c lam rho sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol;  tol; tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;  -1; -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf; inf];
    lb = [-inf*ones(M,1); -inf; tol;  -1; -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf];

end

if strcmp(kernel, 'NEW_bd')
%     % hp = [a1 a2 ... aM c lam rho sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol;  tol; tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;  -1; -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf; inf];
    lb = [-inf*ones(M,1); -inf; tol;  -1; -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf];
end


if strcmp(kernel, 'NEW_new+') || strcmp(kernel, 'NEW_new-')
%     % hp = [a1 a2 ... aM c lam rho sig^2 a0]
%     %lb = [-inf*ones(M,1); -inf; tol;  tol; tol; -inf];
%     lb = [-inf*ones(M,1); -inf; tol;  -1; -inf; -inf];
%     ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf; inf];
    lb = [-inf*ones(M,1); -inf; tol;  -1; -inf];
    ub = [ inf*ones(M,1);  inf; 1-tol; 1; inf];
end
end
