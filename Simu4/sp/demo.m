clear;

% c=parcluster('local');
% c.NumWorkers= 6;
% parpool(6);

M = 3; Lorder = 15; snr  = db2pow(10);

addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);


n = 70;
Nmax = 30000;
Maxrepi = 47;
% noisevar = 0.1;

% for repi = 17%1:Maxrepi
%     data_generation(Nmax, n, M, Lorder, snr, 'expcos', 'gauss',  repi);
% end

N = 500;
kernel = {'DC-dc-sp'};

fprintf('-----------------N = %i-----------------\n',N);
EFIT = []; PFIT = []; COST = []; HP=[]; SIGSQR = [];

for repi = 9
    fprintf('repi = %i: \n',repi);
    d = load([pwd '/Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    data = d.datainfo.data(1:N,:);
    Pi = d.datainfo.Pi(1:N,:);
    Rho = d.datainfo.Rho(1:N,:);
    u = data(:,1);
    y = data(:,2);
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
        EstInfo = rvs_ml_sp(Pi, Rho, y,  n, M, kernel{kk});
        t = toc;
        hyper = EstInfo.hp;
        h0 = EstInfo.h0;
        Oiinvy = EstInfo.Oiinvy;
        O = EstInfo.O;
        ye = EstInfo.yhat;
        cost = EstInfo.cost;
        sigsqr = EstInfo.sigsqr;
        efit = gof(ytrue,ye);
        Ov = CalculateOutputKernelValidation(CalculatePsi(uv, n), CalculatePsi(u, n), M, kernel{kk}, hyper);
        yp = Ov*Oiinvy + h0;
        pfit = gof(yv_true,yp);


        Efit(kk) = efit;
        Sigsqr(kk) = sigsqr;
        Pfit(kk) = pfit;
        Cost(kk) = cost;
        Hp = [Hp; hyper; nan];
        
        fprintf([' finished with times %.4f s\n'],  t);
    end
    EFIT = [EFIT; Efit];
    SIGSQR = [SIGSQR; Sigsqr];
    PFIT = [PFIT; Pfit];
    COST = [COST; Cost];
    HP = [HP Hp];
    fprintf('Completed\n');
end

% p_boxplot(EFIT,10,100,kernel,'EFIT','fit');
% p_boxplot(PFIT,10,100,kernel,'PFIT','fit');
% save('Results/EFIT.mat','EFIT');
% save('Results/PFIT.mat','PFIT');
% save('Results/COST.mat','COST');
% save('Results/HP.mat','HP');