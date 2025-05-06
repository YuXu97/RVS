clear;
% %
c=parcluster('local');
c.NumWorkers= 50;
parpool(50);

M = 3; Lorder = 15; snr  = db2pow(10);

addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);

n = 50;
Nmax = 30000;
Maxrepi = 50;
kernel = {'DC-bd-sp','DC-dc-sp','DC-ob-sp'};%{'DC-dc-sp'};%
N = 500;

Time = []; Time_ini = []; Pfit = []; Cost = []; Sigsqr = [];Hp = []; Hpini = [];
parfor repi = 1:Maxrepi%50%
    warning('off','all')
    %     fprintf('repi = %i: \n',repi);
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    fprintf('repi = %i with N = %i \n',repi,N);
    data = d.datainfo.data(1:N,:);
    Pi = d.datainfo.Pi(1:N,:);
    Rho = d.datainfo.Rho(1:N,:);
    u = data(:,1);
    y = data(:,2);
    Nv_start = N; dN = 5*N;
    datav = d.datainfo.data(Nv_start+1:Nv_start+dN,:);
    uv = datav(:,1); yv_true = d.datainfo.ytrue(Nv_start+n:Nv_start+dN);
    ytrue = d.datainfo.ytrue(n:N);

    time_kernel = zeros(length(kernel),1);
    time_ini_kernel = zeros(length(kernel),1);
    pfit = zeros(length(kernel),1);
    sigsqr = zeros(length(kernel),1);
    cost = zeros(length(kernel),1);
    hp = [];
    hpini = [];
    for kk = 1:length(kernel)
        fprintf(['repi = %i with N = %i and kernel ' kernel{kk}], repi,N);
        tic
        EstInfo = rvs_ml_sp(Pi, Rho, y,  n, M, kernel{kk});
        time_kernel(kk) = toc;
        time_ini_kernel(kk) = EstInfo.time_ini;
        h0 = EstInfo.h0;
        hyper = EstInfo.hp;
        hyperini = EstInfo.hpini;
        Oiinvy = EstInfo.Oiinvy;
        Ov = CalculateOutputKernelValidation(CalculatePsi(uv, n), CalculatePsi(u, n), M, kernel{kk}, hyper);
        yp = Ov*Oiinvy + h0;
        pfit(kk) = gof(yv_true,yp);
        cost(kk) = EstInfo.cost;
        sigsqr(kk) = EstInfo.sigsqr;

        hpini = [hpini; hyperini; nan];
        hp = [hp; hyper; nan];
        fprintf([' finished with times %.4f s\n'],  time_kernel(kk));
    end
    Time = [Time time_kernel];
    Time_ini = [Time_ini time_ini_kernel];
    Pfit = [Pfit pfit];
    Sigsqr = [Sigsqr sigsqr];
    Cost = [Cost cost];
    Hp = [Hp hp];
    Hpini = [Hpini hpini];
end

kns = {'DC-bd-w','DC-decay-w','DC-ob-w'};
nice_boxplot(Pfit',60,100,kns,'','',2);

save([pwd '/Results/Time.mat'], 'Time');
save([pwd '/Results/Time_ini.mat'], 'Time_ini');
save([pwd '/Results/Pfit.mat'], 'Pfit');
save([pwd '/Results/Sigsqr.mat'], 'Sigsqr');
save([pwd '/Results/Cost.mat'], 'Cost');
save([pwd '/Results/Hp.mat'], 'Hp');
save([pwd '/Results/Hpini.mat'], 'Hpini');