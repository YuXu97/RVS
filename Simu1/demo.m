% c=parcluster('local');
% c.NumWorkers= 80;
% parpool(80);
clear;

addpath([pwd '/BasicFunctions']);
addpath([pwd '/conv2fft']);
addpath([pwd '/TuningMethods']);

M = 2; Lorder = 10; snr  = db2pow(5);
n = 80;
Nmax = 5000;
Maxrepi = 80;
% 
% parfor repi = 1:Maxrepi
%     data_generation(Nmax, n, M, Lorder, snr, {'filtered_multisine'}, 'gauss', repi);
% end

N = 400;
fprintf('-----------------N = %i-----------------\n',N);

EFIT = []; PFIT = []; COST = []; HP=[]; SIGSQR = []; OFF = [];

kernel = {'WH-DC'}; method = 'chol';

for repi = 1:Maxrepi
    fprintf('repi = %i: \n',repi);
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    data = d.datainfo.data(1:N,:);
    Nv_start = N; dN = 5*N;
    datav = d.datainfo.data(Nv_start+1:Nv_start+dN,:);
    uv = datav(:,1);
    
    Efit = zeros(1,length(kernel));
    Pfit = zeros(1,length(kernel));
    Cost = zeros(1,length(kernel));
    Sigsqr = zeros(1,length(kernel));
    Off = zeros(1,length(kernel));
    
    Hp = [];
    for kk = 1:length(kernel)
        t = 0;
        fprintf(['repi = %i with kernel ' kernel{kk}], repi)
        if strcmp(kernel{kk}, 'WH-DC')
            nstart = 2*n;
            yv_true = d.datainfo.ytrue(Nv_start+nstart:Nv_start+dN);
            ytrue = d.datainfo.ytrue(nstart:N);
            tic
            EstInfo = rvs_ml(data, n, M, 'WH-DC', method);
            t = toc;
            hyper = EstInfo.hp;
            h0 = EstInfo.h0;
            ye = EstInfo.yhat;
            W  = EstInfo.W;
            cost = EstInfo.cost;
            Psi = EstInfo.Psi;
            sigsqr = EstInfo.sigsqr;
            efit = gof(ytrue,ye);
            Ov = CalculateOutputKernel(CalculatePsi(uv,n), Psi, M, 'WH-DC', hyper, 0);
            yp = Ov*W + h0;
            pfit = gof(yv_true,yp);
            
            Efit(kk) = efit;
            Pfit(kk) = pfit;
            Cost(kk) = cost;
            Sigsqr(kk) = sigsqr;
            Off(kk) = h0;
            Hp = [Hp; hyper; nan];
        elseif strcmp(kernel{kk}, 'DC2-DC')
            nstart = n;
            yv_true = d.datainfo.ytrue(Nv_start+nstart:Nv_start+dN);
            ytrue = d.datainfo.ytrue(nstart:N);
            tic
            EstInfo = rvs_ml_NoHadamard(data, n, M, 'DC2-DC', method);
            t = toc;
            hyper = EstInfo.hp;
            h0 = EstInfo.h0;
            ye = EstInfo.yhat;
            W  = EstInfo.W;
            cost = EstInfo.cost;
            Psi = EstInfo.Psi;
            sigsqr = EstInfo.sigsqr;
            efit = gof(ytrue(n+1:end),ye(n+1:end));
            Phi21f = [];
            Phi21b = [];
            Phi22f = [];
            Phi22b = [];
            Psi1 = CalculatePsi(uv,n);
            Psi2 = Psi;
            for i = 0:n-1
                Phi21b = [Phi21b Psi2(:,1:n-i)];
                if i ~= 0
                    Phi22b = [Phi22b 2*Psi2(:,i+1:n)];
                else
                    Phi22b = [Phi22b Psi2(:,i+1:n)];
                end
                
                Phi21f = [Phi21f Psi1(:,1:n-i)];
                if i ~= 0
                    Phi22f = [Phi22f 2*Psi1(:,i+1:n)];
                else
                    Phi22f = [Phi22f Psi1(:,i+1:n)];
                end
            end
            
            Phi2b = Phi21b.*Phi22b;
            Phi2f = Phi21f.*Phi22f;
            
            [TI,TJ] = meshgrid(1:n, 1:n);
            indxi = []; indxj = [];
            for i = 0:n-1
                indxi = [indxi; diag(TI,i)];
            end
            for i = 0:n-1
                indxj = [indxj; diag(TJ,i)];
            end
            indx = [indxj,indxi];
            indx = indx*[cos(pi/4), -sin(pi/4);sin(pi/4), cos(pi/4)];
            
            Ov = CalculateOutputKernel_NoHadamard(Psi1, Psi2, Phi2f, Phi2b, indx, M, 'DC2-DC', hyper, 0);
            yp = Ov*W + h0;
            pfit = gof(yv_true(n+1:end),yp(n+1:end));
            
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
    
    %WHDC
save('Results/EFIT.mat','EFIT');
save('Results/PFIT.mat','PFIT');
save('Results/COST.mat','COST');
save('Results/SIGSQR.mat','SIGSQR');
save('Results/OFF.mat','OFF');
save('Results/HP.mat','HP');
end


p_boxplot(EFIT,20,100,kernel,'EFIT','fit');
p_boxplot(PFIT,10,100,kernel,'PFIT','fit');
