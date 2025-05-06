clear;
% 
c=parcluster('local');
c.NumWorkers= 80;
parpool(80);

addpath([pwd '/BasicFunctions']);
addpath([pwd '/conv2fft']);
addpath([pwd '/TuningMethods']);

M = 2; Lorder = 30; snr  = db2pow(10); 
n = 80;
Nmax = 8000;

% data_generation_accuracy_N(Nmax, n, M, Lorder, snr, {'white',[0,1]}, 'gauss', 0);
d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(0) '.mat']);
Nrange = 2200:200:2600;

EFIT = []; PFIT = [];

kernel = 'DC_wh_bd';

method = 'chol';

for repi = 1:length(Nrange)
    N = Nrange(repi);
    fprintf('-----------------N = %i-----------------\n',N);
    data = d.datainfo.data(1:N,:);
    Nv_start = N; dN = N;
    datav = d.datainfo.data(Nv_start+1:Nv_start+dN,:);
    uv = datav(:,1);
    nstart = 2*n;
    yv_true = d.datainfo.ytrue(Nv_start+nstart:Nv_start+dN);
    ytrue = d.datainfo.ytrue(nstart:N);
    tic;
    EstInfo = rvs_ml(data, n, M, kernel, method);
    t = toc;
    hyper = EstInfo.hp;
    h0 = EstInfo.h0;
    ye = EstInfo.yhat;
    W  = EstInfo.W;
    cost = EstInfo.cost;
    Psi = EstInfo.Psi;
    sigsqr = EstInfo.sigsqr;
    efit = gof(ytrue,ye);
    Ov = CalculateOutputKernel(CalculatePsi(uv,n), Psi, M, kernel, hyper, 0);
    yp = Ov*W + h0;
    pfit = gof(yv_true,yp);
    fprintf([' finished with times %.4f s\n'],  t);
    EFIT = [EFIT efit];
    PFIT = [PFIT pfit];
end


save('Results/PFIT_2600.mat', 'PFIT');
save('Results/EFIT_2600.mat', 'EFIT');
