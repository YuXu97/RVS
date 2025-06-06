clear;
addpath([pwd '/Regular'])
addpath([pwd '/SP'])
warning('off','all')
d = load('data_N30000_repi=49');

kernel = {'DC-bd-sp'};

N = 400;
M = 3;
n = 50;
data = d.datainfo.data(1:N,:);
u = data(:,1);
y = data(:,2);


 hyper =   [8;-0.888888888888889;0.888888888888889;0.881111111111111];
%  [3.67816200449581;
%  4.38351123173726;
%  -7.57597181481623;
%  0.990446431503013];
%   hyper =     [3.678162004495810;4.383511231737260;-7.575971814816230;0.990446431503013];
%  hyper =  [3.844490851526519;4.541382020729821;-7.459890162867881;0.906299079415358];
%  hyper = [1;1;1;0.99];

tol = sqrt(eps);
lb = [-inf*ones(M,1); tol];
ub = [inf*ones(M,1); 1-tol];

%regular
Psi = CalculatePsi(u, n);
ycut = data(n:N,2);
method = 'chol';
[obj_reg, qdr_reg, lgdt_reg] = nglglklhd_reg(hyper, Psi, ycut, kernel{1}, M, method);
% ff_reg = @(x)nglglklhd_reg(x, Psi, ycut, kernel{1}, M, method);
% % options = optimoptions('fmincon', 'Display', 'off',  'Algorithm', 'active-set');
% % hp_reg = fmincon(ff_reg,hyper,[],[],[],[],lb,ub,[],options);
% options = optimset('Display', 'off');
% hp_reg = fminsearch(ff_reg,hyper,options);
% 
% obj_reg_final = nglglklhd_reg(hp_reg, Psi, ycut, kernel{1}, M, method);



%sp
Pi = d.datainfo.Pi(n:N,:);
Rho = d.datainfo.Rho(1:n,:);
ycut = data(n:N,2);

[~,r] = size(Pi);
pb = r;
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

[obj_sp, qdr_sp, lgdt_sp] = nglglklhd_sp(hyper, Pi, Rho, ycut, kernel{1}, M, n, CM, Cf, indcum);
% ff_sp = @(x)nglglklhd_sp(x,  Pi, Rho, ycut, kernel{1}, M, n, CM, Cf, indcum);
% % options = optimoptions('fmincon', 'Display', 'off',  'Algorithm', 'active-set');
% % hp_sp = fmincon(ff_sp,hyper,[],[],[],[],lb,ub,[],options);
% 
% options = optimset('Display', 'off');
% hp_sp = fminsearch(ff_sp,hyper,options);
% 
% obj_sp_final = nglglklhd_sp(hp_sp, Pi, Rho, ycut, kernel{1}, M, n, CM, Cf, indcum);

obj_sp - obj_reg