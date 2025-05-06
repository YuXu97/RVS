clear;
% %
% c=parcluster('local');
% c.NumWorkers= 6;
% parpool(6);

M = 3; Lorder = 15; snr  = db2pow(10);

addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);
addpath([pwd '/EGRSS_m']);

n = 50;
Nmax = 30000;
Maxrepi = 50;
kernel = {'DC-bd-sp','DC-dc-sp','DC-ob-sp'};
Nrange = 1000:1000:8000;

warning('off','all')
for repi = 1:Maxrepi
    %     fprintf('repi = %i: \n',repi);
    Time = []; Time_ini = []; Pfit = []; Cost = []; Sigsqr = [];
    Hp = [];
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    for iterN = 1:length(Nrange)    
        N = Nrange(iterN);
        %fprintf('repi = %i with N = %i \n',repi,N);
        data = d.datainfo.data(1:N,:);
        Pi = d.datainfo.Pi(1:N,:);
        Rho = d.datainfo.Rho(1:N,:);
        u = data(:,1);
        y = data(:,2);
        Nv_start = N; dN = 8000;%N;
        datav = d.datainfo.data(Nv_start+1:Nv_start+dN,:);
        uv = datav(:,1); yv_true = d.datainfo.ytrue(Nv_start+n:Nv_start+dN);
        ytrue = d.datainfo.ytrue(n:N);

        time_kernel = zeros(length(kernel),1);
        time_ini_kernel = zeros(length(kernel),1);
        pfit = zeros(length(kernel),1);
        sigsqr = zeros(length(kernel),1);
        cost = zeros(length(kernel),1);
        hp = [];
        for kk = 1:length(kernel)
            fprintf(['repi = %i with N = %i and kernel ' kernel{kk}], repi,N);
            tic
            EstInfo = rvs_ml_smsp(Pi, Rho, y,  n, M, kernel{kk});
            time_kernel(kk) = toc;
            time_ini_kernel(kk) = EstInfo.time_ini;
            h0 = EstInfo.h0;
            hyper = EstInfo.hp;
            Oiinvy = EstInfo.Oiinvy;
            Ov = CalculateOutputKernelValidation(CalculatePsi_ir(uv, n), CalculatePsi_ir(u, n), M, kernel{kk}, hyper);
            yp = Ov*Oiinvy + h0;
            pfit(kk) = gof(yv_true,yp);
            cost(kk) = EstInfo.cost;
            sigsqr(kk) = EstInfo.sigsqr;

            hp = [hp; hyper; nan];
            fprintf([' finished with times %.4f s\n'],  time_kernel(kk));
        end
        Time = [Time time_kernel];
        Time_ini = [Time_ini time_ini_kernel];
        Pfit = [Pfit pfit];
        Sigsqr = [Sigsqr sigsqr];
        Cost = [Cost cost];
        Hp = [Hp hp];

    end

    save([pwd '/Results/time/Time' int2str(repi) '.mat'], 'Time');
    save([pwd '/Results/time_ini/Time_ini' int2str(repi) '.mat'], 'Time_ini');
    save([pwd '/Results/Pfit/Pfit' int2str(repi) '.mat'], 'Pfit');
    save([pwd '/Results/Sigsqr/Sigsqr' int2str(repi) '.mat'], 'Sigsqr');
    save([pwd '/Results/Cost/Cost' int2str(repi) '.mat'], 'Cost');
    save([pwd '/Results/Hp/Hp' int2str(repi) '.mat'], 'Hp');
end

Time = zeros(length(kernel),length(Nrange));
Time_ini = zeros(length(kernel),length(Nrange));
for repi = 1:Maxrepi
    t = load([pwd '/Results/time/Time' int2str(repi) '.mat']);
    t_ini = load([pwd '/Results/time_ini/Time_ini' int2str(repi) '.mat']);
    time = t.Time;
    time_ini = t_ini.Time_ini;
    Time = Time + time;
    Time_ini = Time_ini + time_ini;
end
Time = Time/Maxrepi;
Time_ini = Time_ini/Maxrepi;
figure(1)
nice_plot(Nrange,Time,2); 
kns = {'DC-bd-w','DC-decay-w','DC-ob-w'};
legend(kns)
figure(2)
nice_plot(Nrange,Time_ini,2);
legend(kns)