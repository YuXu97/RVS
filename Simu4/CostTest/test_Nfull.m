clear;
addpath([pwd '/Regular_N'])
addpath([pwd '/SMSP'])

warning('off','all')
d = load('data_N30000_repi=50.mat');

kernel = {'DC-ob-sp'};

N = 400;
M = 3;
n = 50;
data = d.datainfo.data(1:N,:);
u = data(:,1);
y = data(:,2);


hyper =  [8;-0.888888888888889;0.888888888888889;0.881111111111111];%[1;1;1;0.9];%

tol = sqrt(eps);
lb = [-inf*ones(M,1); tol];
ub = [inf*ones(M,1); 1-tol];

%regular_N
Psi = CalculatePsi_ir(u, n);
ycut = data(n:N,2);
method = 'chol';
[obj_regN, qdr_regN, lgdt_regN] = nglglklhd_regN(hyper, Psi, ycut, kernel{1}, M, method);
% ff_regN = @(x)nglglklhd_regN(x, Psi, ycut, kernel{1}, M, method);
% options = optimset('Display', 'off');
% hp_regN = fminsearch(ff_regN,hyper,options);
% 
% obj_regN_final = nglglklhd_regN(hp_regN, Psi, ycut, kernel{1}, M, method);

%smsp
Pi = d.datainfo.Pi(1:N,:);
Rho = d.datainfo.Rho(1:N,:);
Psi_ir = CalculatePsi_ir(u, n);
ycut = data(n:N,2);

p = 1;
[~,r] = size(Pi);
pb = p + r;

indm = zeros(M,1);
for m = 1:M
    indm(m) = factorial(pb+m-1)/factorial(m)/factorial(pb-1);
end
indcum = cumsum(indm); indcum = [0;indcum];
CM = zeros(sum(indm),pb);
Cf = zeros(sum(indm),1);
for m = 1:M
    C = MultiPermuations(m,pb);
    CM(indcum(m)+1:indcum(m+1),:) = C;
    for i = 1:indm(m)
        Cf(indcum(m)+i) = multinomial(m,C(i,:));
    end
end

[obj_smsp, qdr_smsp, lgdt_smsp] = nglglklhd_smsp(hyper, Pi, Rho, Psi_ir, ycut, kernel{1}, M, n, CM, Cf, indcum);

% ff_smsp = @(x)nglglklhd_smsp(x, Pi, Rho, Psi_ir, ycut, kernel{1}, M, n, CM, Cf, indcum);
% options = optimset('Display', 'off');
% hp_smsp = fminsearch(ff_smsp,hyper,options);
% 
% obj_smsp_final = nglglklhd_smsp(hp_smsp, Pi, Rho, Psi_ir, ycut, kernel{1}, M, n, CM, Cf, indcum);


obj_regN - obj_smsp
