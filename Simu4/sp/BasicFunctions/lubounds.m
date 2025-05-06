function [lb,ub] = lubounds(kernel, M)
tol = sqrt(eps);

if strcmp(kernel, 'TC-bd') || strcmp(kernel, 'TC-tc') ...
        || strcmp(kernel, 'DC-bd-sp') 
     % hp = [c1 ... CM lam]
    lb = [-inf*ones(M,1);  tol];
    ub = [inf*ones(M,1);  1-tol];

end

if strcmp(kernel, 'DC-dc-sp') || strcmp(kernel, 'DC-ob-sp')
    lb = [tol; -inf*ones(M-1,1);  tol];
    ub = [inf; inf*ones(M-1,1);  1-tol];
end

if strcmp(kernel, 'DC-bd')
    lb = [-inf*ones(M,1);  tol; tol];
    ub = [inf*ones(M,1);  1-tol; tol];
end

if  strcmp(kernel, 'DC-dc') || strcmp(kernel, 'DC-ob')
    % hp = [c1 ... CM lam]
    lb = [tol; -inf*ones(M-1,1);  tol; tol];
    ub = [inf; inf*ones(M-1,1);  1-tol; 1];
end

end
