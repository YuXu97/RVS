%  c=parcluster('local');
%  c.NumWorkers= 85;
%  parpool(85);
% delete(gcp('nocreate'));
clear;
addpath([pwd '/BasicFunctions']);
addpath([pwd '/conv2fft']);
addpath([pwd '/TuningMethods']);

M = 3; Lorder = 30; snr  = db2pow(5);  
n = 80;
Nmax = 4000;
Maxrepi = 200;
%  
% parfor i=1:length(lack_list)%repi = 1:Maxrepi%
%              repi = lack_list(i);
%     data_generation(Nmax, n, M, Lorder, snr, {'white',[0,1]}, 'gauss', repi);
% end

N = 400;
fprintf('-----------------N = %i-----------------\n',N);

EFIT = []; PFIT = []; COST = []; HP=[]; SIGSQR = []; OFF = [];

% kernel = {'DC_mpoly','DC_wh_bd','DC_wh_dc','DC_wh_ob','DC_wh_oba'};
kernel = {'DC_mpoly','DC_wh_bd','DC_wh_dc','DC_wh_oba'};
% kernel = {'DC_wh_bd','DC_wh_oba'};

method = 'chol';

parfor repi = 1:Maxrepi 
    fprintf('repi = %i: \n',repi);
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    data = d.datainfo.data(1:N,:);
    Nv_start = N; dN = 5*N;
    datav = d.datainfo.data(Nv_start+1:Nv_start+dN,:);
    uv = datav(:,1);
    yv_true = d.datainfo.ytrue(Nv_start+n:Nv_start+dN);
    ytrue = d.datainfo.ytrue(n:N);
    
    Efit = zeros(1,length(kernel));
    Pfit = zeros(1,length(kernel));
    Cost = zeros(1,length(kernel));
    Sigsqr = zeros(1,length(kernel));
    Off = zeros(1,length(kernel));
    
    Hp = [];
    for kk = 1:length(kernel)
        fprintf(['repi = %i with kernel ' kernel{kk}], repi)
        
        if strcmp(kernel{kk}, 'DC_mpoly')
            nstart = n;
            yv_true = d.datainfo.ytrue(Nv_start+nstart:Nv_start+dN);
            ytrue = d.datainfo.ytrue(nstart:N);
            tic
            EstInfo = rvs_ml(data, n, M, kernel{kk}, method);
            t = toc;
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
            %pfit = gof(yv_true,yp);
            pfit = gof(yv_true(n+1:end),yp(n+1:end));
            
            Efit(kk) = efit;
            Pfit(kk) = pfit;
            Cost(kk) = cost;
            Sigsqr(kk) = sigsqr;
            Off(kk) = h0;
            Hp = [Hp; hyper; nan];
        else
            nstart = 2*n;
            yv_true = d.datainfo.ytrue(Nv_start+nstart:Nv_start+dN);
            ytrue = d.datainfo.ytrue(nstart:N);
            tic
            EstInfo = rvs_ml(data, n, M, kernel{kk}, method);
            t = toc;
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
            Pfit(kk) = pfit;
            Cost(kk) = cost;
            Sigsqr(kk) = sigsqr;
            Off(kk) = h0;
            Hp = [Hp; hyper; nan];
        end
        fprintf([' finished with times %.4f s\n'],  t);
    end
    EFIT = [EFIT; Efit];
    PFIT = [PFIT; Pfit];
    COST = [COST; Cost];
    SIGSQR = [SIGSQR; Sigsqr];
    OFF = [OFF; Off];
    HP = [HP Hp];
    fprintf('Completed\n');
end

% nice_boxplot(EFIT,50,100,kernel,'EFIT','fit',2);
nice_boxplot(PFIT,20,100,kernel,'PFIT','fit',2);
save('Results/EFIT.mat','EFIT');
save('Results/PFIT.mat','PFIT');
save('Results/COST.mat','COST');
save('Results/SIGSQR.mat','SIGSQR');
save('Results/OFF.mat','OFF');
save('Results/HP.mat','HP');