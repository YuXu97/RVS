clear;
% 
% c=parcluster('local');
% c.NumWorkers= 80;
% parpool(80);

addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);

M = 2; n = 20;
Nmax = 1000;
Maxrepi = 80;
% 
% for repi = 1:Maxrepi
%     data_generation(Nmax, n, {'white',[0,1]}, 0.95, 'gauss', repi);
% end

Efits_poly = [];Efits_mpoly = [];Efits_LNL= [];
Costs_poly = [];Costs_mpoly = [];Costs_LNL = [];
N = 100;


    
fprintf('-----------------N = %i-----------------\n',N);
Efit_poly = [];Efit_mpoly = [];Efit_LNL= [];
Pfit_poly = [];Pfit_mpoly = [];Pfit_LNL= [];
Cost_poly = [];Cost_mpoly = [];Cost_LNL = [];

for repi = 48:48

d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
data = d.datainfo.data(1:N,:);
datav = d.datainfo.data(N+1:2*N,:);
uv = datav(:,1); yv_true = d.datainfo.ytrue(N+n:2*N);

%--------------------------------training----------------------------------%
kernel = {'DC_poly','DC_mpoly','DC_LNL'}; method = 'svd';

fprintf('repi = %i: \n',repi);
tic;
EstInfo_poly = rvs_ml(data, n, M, kernel{1}, method);
EstInfo_mpoly = rvs_ml(data, n, M, kernel{2}, method);
% EstInfo_LNL = rvs_ml(data, n, M, kernel{3}, method);
toc;
fprintf('Completed\n');

ytrue = d.datainfo.ytrue(n:N);

hyper_poly = EstInfo_poly.hp;
ye_poly = EstInfo_poly.yhat;
W_poly  = EstInfo_poly.W;
cost_poly = EstInfo_poly.cost; 
Psi_poly = EstInfo_poly.Psi;
efit_poly = gof(ytrue,ye_poly);
Efit_poly = [Efit_poly;efit_poly];
Cost_poly = [Cost_poly;cost_poly];
Ov_poly = CalculateOutputKernel(CalculatePsi(uv,n),Psi_poly, M, kernel{1}, hyper_poly, 0);
yp_poly = Ov_poly*W_poly + hyper_poly(end);
pfit_poly = gof(yv_true,yp_poly);
Pfit_poly = [Pfit_poly;pfit_poly];


hyper_mpoly = EstInfo_mpoly.hp;
ye_mpoly = EstInfo_mpoly.yhat;
W_mpoly  = EstInfo_mpoly.W;
cost_mpoly = EstInfo_mpoly.cost; 
Psi_mpoly = EstInfo_mpoly.Psi;
efit_mpoly = gof(ytrue,ye_mpoly);
Efit_mpoly = [Efit_mpoly;efit_mpoly];
Cost_mpoly = [Cost_mpoly;cost_mpoly];
Ov_mpoly = CalculateOutputKernel(CalculatePsi(uv,n),Psi_mpoly, M, kernel{2}, hyper_mpoly, 0);
yp_mpoly = Ov_mpoly*W_mpoly + hyper_mpoly(end);
pfit_mpoly = gof(yv_true,yp_mpoly);
Pfit_mpoly = [Pfit_mpoly;pfit_mpoly];

% hyper_LNL = EstInfo_LNL.hp;
% ye_LNL = EstInfo_LNL.yhat;
% W_LNL  = EstInfo_LNL.W;
% cost_LNL = EstInfo_LNL.cost;
% Psi_LNL = EstInfo_LNL.Psi;
% efit_LNL = gof(ytrue,ye_LNL);
% Efit_LNL = [Efit_LNL;efit_LNL];
% Cost_LNL= [Cost_LNL;cost_LNL];
% Ov_LNL = CalculateOutputKernel(CalculatePsi(uv,n),Psi_LNL, M, kernel{3}, hyper_LNL, 0);
% yp_LNL = Ov_LNL*W_LNL + hyper_LNL(end);
% pfit_LNL = gof(yv_true,yp_LNL);
% Pfit_LNL = [Pfit_LNL;pfit_LNL];
end

% save('Efit_poly.mat', 'Efit_poly');
% save('Pfit_poly.mat', 'Pfit_poly');
% save('Cost_poly.mat', 'Cost_poly');
% 
% save('Efit_mpoly.mat', 'Efit_mpoly');
% save('Pfit_mpoly.mat', 'Pfit_mpoly');
% save('Cost_mpoly.mat', 'Cost_mpoly');

% save('Efit_LNL.mat', 'Efit_LNL');
% save('Pfit_LNL.mat', 'Pfit_LNL');
% save('Cost_LNL.mat', 'Cost_LNL');
% 
% p_boxplot([Efit_poly Efit_mpoly],0,100,{'poly','mpoly'},'','');
% 
% p_boxplot([Pfit_poly Pfit_mpoly],0,100,{'poly','mpoly'},'','');

% p_boxplot([Efit_LNL Efit_mpoly],0,100,{'LNL','mpoly'},'','');
% 
% p_boxplot([Pfit_LNL Pfit_mpoly],0,100,{'LNL','mpoly'},'','');