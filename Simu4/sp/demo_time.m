clear;
% %
% c=parcluster('local');
% c.NumWorkers= 6;
% parpool(6);

M = 3; Lorder = 15; snr  = db2pow(10);

addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);


n = 40;
Nmax = 60000;
Maxrepi = 1;
noisevar = 0.1;
% 
% for repi = 1:Maxrepi
%     data_generation(Nmax, n, M, Lorder, snr, 'expcos', 'gauss', repi);
% end

% N = 500;

kernel = {'DC-bd-sp','DC-dc-sp','DC-ob-sp'};

EFIT = []; PFIT = []; COST = []; HP=[]; SIGSQR = [];TIME = [];
Nrange = 1000:1000:10000; time_rp = 1;
time = zeros(length(Nrange),1);
for iter = 1:length(Nrange)
    N = Nrange(iter);
    fprintf('-----------------N = %i-----------------\n',N);
    repi = 1;
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    data = d.datainfo.data(1:N,:);
    Pi = d.datainfo.Pi(1:N,:);
    Rho = d.datainfo.Rho(1:N,:);
    u = data(:,1);
    y = data(:,2);
    Nv_start = N; dN = N;
    datav = d.datainfo.data(Nv_start+1:Nv_start+dN,:);
    uv = datav(:,1); yv_true = d.datainfo.ytrue(Nv_start+n:Nv_start+dN);
    ytrue = d.datainfo.ytrue(n:N);
    Efit = zeros(1,length(kernel));
    Sigsqr = zeros(1,length(kernel));
    Pfit = zeros(1,length(kernel));
    Cost = zeros(1,length(kernel));
    Hp = [];
    Time = zeros(1,length(kernel));
    for kk = 1:length(kernel)
        fprintf(['N = %i with kernel ' kernel{kk}], N)
        t = zeros(time_rp,1);
        for rp = 1:time_rp
            tic
            EstInfo = rvs_ml_sp(Pi, Rho, y,  n, M, kernel{kk});
            t(rp) = toc;
        end

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

        Time(kk) = mean(t);
        Efit(kk) = efit;
        Sigsqr(kk) = sigsqr;
        Pfit(kk) = pfit;
        Cost(kk) = cost;
        Hp = [Hp; hyper; nan];

        fprintf(' finished with average time %.4f s\n',  Time(kk));
    end
    EFIT = [EFIT; Efit];
    SIGSQR = [SIGSQR; Sigsqr];
    PFIT = [PFIT; Pfit];
    COST = [COST; Cost];
    HP = [HP Hp];
    TIME = [TIME; Time]
    fprintf('Completed\n');
end

% p_boxplot(EFIT,10,100,kernel,'EFIT','fit');
% p_boxplot(PFIT,10,100,kernel,'PFIT','fit');
figure(3)
plot(Nrange,TIME); grid on; xlabel('N'); ylabel('avtime');
legend(kernel)
% plot(Nrange,time); hold on;
% plot(Nrange,timedc); xlabel('N'); ylabel('avtime');legend('DC-bd','DC-dc'); grid on;
save('Results/EFIT.mat','EFIT');
save('Results/PFIT.mat','PFIT');
save('Results/COST.mat','COST');
save('Results/HP.mat','HP');
save('Results/time.mat','TIME');