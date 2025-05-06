clear;
% c=parcluster('local');
% c.NumWorkers= 5;
% parpool(5);

addpath('generatedata/');
addpath(genpath('No_Identifiability_gmle'))
N= 400;
n = 80;
par.finalTime       = N;               % Data length
par.resampling      = 3;               % 1 - multinom, 2 - systematic, 3 - stratified
par.Np              = 15;              % Number of particles
par.numMCMC         = 10000;           % Number of MCMC iterations
par.burnin = 5000;                     % MCMC burnin. Also, stops rescaling after this many iterations.
par.input           = 1;               % Use this to indicate whether or not we should use an input signal when generating data
par.scaling         = 'h2';            % How to resolve the scale ambiguity, either 'range', 'H2' or 'none'

%
par.dynprior        = 'ard';
par.ardnx           = 5;

Maxrepi = 1;

% for repi = 1:Maxrepi
%     testsystem_wienergp_2d(N, par, repi);
% end

PFIT_vs = [];
PFIT_nlwh = [];

for repi = 1:Maxrepi
d = load(['databank/' 'data_' int2str(repi) '.mat']);
data = [d.data.u(1:N)',d.data.y(1:N)'];

data = flip(data,2);
order = [2 2 1];
pl =idPolynomial1D('Degree',2);
sys = nlhw(data, order, [], pl);
y_nlhw = sim(sys,d.data.uv');
PFIT_nlwh = [PFIT_nlwh gof(d.data.yvtrue, y_nlhw)];
% 
M = 2;
EstInfo = rvs_ml([d.data.u(1:N)',d.data.y(1:N)'], n, M, 'DC-bd', 'chol');
hyper = EstInfo.hp;
h0 = EstInfo.h0;
W  = EstInfo.W;
Psi = EstInfo.Psi;
Ov = CalculateOutputKernel(CalculatePsi(d.data.uv(1:N)',n), Psi, M, 'DC-bd', hyper, 0);
yp = Ov*W + h0;
PFIT_vs = [PFIT_vs gof(d.data.ytrue(N+n:2*N),yp)];
end

PFITs = [PFIT_nlwh', PFIT_vs'];

save('PFITs.mat','PFITs')


