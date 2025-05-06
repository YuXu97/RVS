clear;
%
% c=parcluster('local');
% c.NumWorkers= 75;
% parpool(75);


addpath('generatedata/')
addpath('pmcmc/')
addpath('helpers/')
addpath(genpath('gpml-matlab-v3.1-2010-09-27'))
addpath(genpath('No_Identifiability_gmle'))

N= 400;

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

% parfor repi = 1:160
%     testsystem_wienergp_2d(N, par, repi);
% end


M = 2;
n = 70;
EFIT = [];
PFIT = [];
for repi = 1:1
    fprintf('repi = %i: \n',repi);
    d = load(['databank/' 'data_' int2str(repi) '.mat']);
    tic;
    EstInfo = rvs_ml([d.data.u(1:N)',d.data.y(1:N)'], n, M, 'DC-bd', 'chol');
    t = toc;
    hyper = EstInfo.hp;
    h0 = EstInfo.h0;
    W  = EstInfo.W;
    ye = EstInfo.yhat;
    cost = EstInfo.cost;
    sigsqr = EstInfo.sigsqr;


    ytrue_e = d.data.y'-d.data.e(1:N)';
    efit = gof(ytrue_e(n:N),ye);
    EFIT = [EFIT; efit];
    Psi = EstInfo.Psi;
    Ov = CalculateOutputKernel(CalculatePsi(d.data.uv(1:N)',n), Psi, M, 'DC-bd', hyper, 0);
    yp = Ov*W + h0;
    ytrue_v = d.data.yv'-d.data.e(N+1:2*N)';
    pfit = gof(ytrue_v(n:N),yp)
    PFIT = [PFIT; pfit];
    fprintf([' finished with time %.4f s\n'],  t);
end
save('EFIT_volterra.mat','EFIT');
save('PFIT_volterra.mat','PFIT');

% figure(1)
% plot(ytrue_e(n:N)); hold on;
% plot(ye);grid on;
% legend('TRUE','ESTIMATE')
% title('efit')
%
% figure(2)
% plot(ytrue_v(n:N)); hold on;
% plot(yp);grid on;
% legend('TRUE','ESTIMATE')
% title('pfit')

