function [lb,ub] = lubounds(kernel, M)
tol = sqrt(eps);

if strcmp(kernel, 'WH-DC')
    % hp = [c1 c2 lam1 lam2 lam3 rho1 rho2 rho3]

    lb = [-inf*ones(2,1);  tol; tol; tol; tol; tol; tol];
    ub = [inf*ones(2,1); 1-tol; 1-tol; 1-tol;  1;  1; 1];

end

if  strcmp(kernel, 'DC2-DC')
    % hp = [c1 c2 lam1 lam2 lam3 rho1 rho2 rho3]

    lb = [-inf*ones(2,1);  tol; tol; tol; tol; tol; tol];
    ub = [inf*ones(2,1); 1-tol; 1-tol; 1-tol;  1;  1; 1];

end
end
