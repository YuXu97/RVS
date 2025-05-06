% c=parcluster('local');
% c.NumWorkers= 50;
% parpool(50);
clear;
% 
% c=parcluster('local');
% c.NumWorkers= 80;
% parpool(80);


M = 2; Lorder = 15; snr  = db2pow(1000);

% addpath([pwd '/BasicFunctions']);
% addpath([pwd '/TuningMethods']);

n = 70;
Nmax = 4000;
Maxrepi = 1;
noisevar = 0.1;

% for repi = 1:Maxrepi
%     data_generation(Nmax, n, M, Lorder, snr, {'white',[0,1]}, 0.95, 'gauss', noisevar, repi);
% end

N = 400;
kernel = {'DC-bd'};
% kernel = {'DC-ob','DC-eig'};
% kernel = {'DC-bd','DC-ob','DC-ob-ex'};
method = 'chol';

fprintf('-----------------N = %i-----------------\n',N);
EFIT = []; PFIT = []; COST = []; HP=[]; SIGSQR = [];
for repi = 1:Maxrepi
    fprintf('repi = %i: \n',repi);
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    data = d.datainfo.data(1:N,:);
    Nv_start = N; dN = 5*N;
    datav = d.datainfo.data(Nv_start+1:Nv_start+dN,:);
    uv = datav(:,1); yv_true = d.datainfo.ytrue(Nv_start+n:Nv_start+dN);
    ytrue = d.datainfo.ytrue(n:N);
    Efit = zeros(1,length(kernel));
    Sigsqr = zeros(1,length(kernel));
    Pfit = zeros(1,length(kernel));
    Cost = zeros(1,length(kernel));
    Hp = [];
    for kk = 1:length(kernel)
        fprintf(['repi = %i with kernel ' kernel{kk}], repi)
        tic
        EstInfo = rvs_ml(data, n, M, kernel{kk}, method);
        hyper = EstInfo.hp;
        h0 = EstInfo.h0;
        ye = EstInfo.yhat;
        W  = EstInfo.W;
        cost = EstInfo.cost;
        Psi = EstInfo.Psi;
        sigsqr = EstInfo.sigsqr;
        efit = gof(ytrue,ye);

        Ov = CalculateOutputKernel(CalculatePsi(uv,n), Psi, M, kernel{kk}, hyper, 0);
        yp = Ov*W + h0;
        pfit = gof(yv_true,yp);

        Efit(kk) = efit;
        Sigsqr(kk) = sigsqr;
        Pfit(kk) = pfit;
        Cost(kk) = cost;
        Hp = [Hp; hyper; nan];
        t = toc;
        fprintf([' finished with time %.4f s\n'],  t);
    end
    EFIT = [EFIT; Efit];
    SIGSQR = [SIGSQR; Sigsqr];
    PFIT = [PFIT; Pfit];
    COST = [COST; Cost];
    HP = [HP Hp];
    fprintf('Completed\n');
end

p_boxplot(EFIT,20,100,kernel,'EFIT','fit');
p_boxplot(PFIT,30,100,'','PFIT','fit');
save('Results/EFIT.mat','EFIT');
save('Results/PFIT.mat','PFIT');
save('Results/COST.mat','COST');
save('Results/HP.mat','HP');