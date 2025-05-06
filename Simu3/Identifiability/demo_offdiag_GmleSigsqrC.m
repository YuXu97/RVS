

% c=parcluster('local');
% c.NumWorkers= 80;
% parpool(80);
clear;


M = 2; Lorder = 15; snr  = 1;

addpath([pwd '/BasicFunctions']);
addpath([pwd '/TuningMethods']);

n = 20;
Nmax = 3000;
Maxrepi = 80;
poles = zeros(Maxrepi,1);
% for repi = 1:Maxrepi
%     data_generation(Nmax, n, M, Lorder, snr, {'white',[0,1]}, 0.95, 'gauss', repi);
% end


N = 80;

fprintf('-----------------N = %i-----------------\n',N);

Hp = [];Poly = [];

EFIT = [];
PFIT = [];
GFIT = [];
PLFIT = [];
COST = [];
SIGSQR = [];
% kernel = {'DC_eig'};
kernel = {'DC_bd','DC_opt','DC_dc'};
%kernel = {'NEW_bd','NEW_new'};
% ID = load('Results/CHK_M2/ind.mat');
% ind = ID.ind;

parfor repi =1:Maxrepi % 1:20 %
    
    d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    data = d.datainfo.data(1:N,:);
    Psi = CalculatePsi(data(:,1), n);
    ytrue = d.datainfo.ytrue(n:N);
    gtrue = d.datainfo.LinearSystemImpulseResponse(1:n);
    poly_true = d.datainfo.PolyCoefficients;
    Nv_start = N; dN = 5*N;
    datav = d.datainfo.data(Nv_start+1:Nv_start+dN,:);
    uv = datav(:,1); yv_true = d.datainfo.ytrue(Nv_start+n:Nv_start+dN);
    
    %--------------------------------training----------------------------------%
    
    
    method = 'svd';
    
    fprintf('repi = %i: \n',repi);
    Efit = [];
    Pfit = [];
    Gfit = [];
    Plfit = [];
    Cost = [];
    Sigsqr = [];
    tic;
    for kk = 1:length(kernel)
        
        if strcmp(kernel{kk}, 'TC_opt') || strcmp(kernel{kk}, 'TC_tc')...
                || strcmp(kernel{kk}, 'DC_opt') || strcmp(kernel{kk}, 'DC_dc') ...
                || strcmp(kernel{kk}, 'NEW_new') || strcmp(kernel{kk}, 'DC_eig')
            EstInfo1 = rvs_ml_GmleSigsqrC(data,  n,  M, [kernel{kk},'+'], method);
            EstInfo2 = rvs_ml_GmleSigsqrC(data,  n,  M, [kernel{kk},'-'], method);
            if EstInfo1.cost <  EstInfo2.cost
                EstInfo = EstInfo1;
                hyper = EstInfo.hp;
                Ov = CalculateOutputKernel(CalculatePsi(uv,n),Psi, M, [kernel{kk},'+'], hyper, 0);
            else
                EstInfo = EstInfo2;
                hyper = EstInfo.hp;
                Ov = CalculateOutputKernel(CalculatePsi(uv,n),Psi, M, [kernel{kk},'-'], hyper, 0);
            end
        else
            EstInfo = rvs_ml_GmleSigsqrC(data,  n,  M, kernel{kk}, method);
            hyper = EstInfo.hp;
            Ov = CalculateOutputKernel(CalculatePsi(uv,n),Psi, M, kernel{kk}, hyper, 0);
        end
        
        
        Ov = exp(hyper(end))*Ov;
        
        poly = hyper(1:M);
        
        ghat = EstInfo.ghat;
        b = EstInfo.b;
        ye = EstInfo.yhat;
        W  = EstInfo.W;
        cost = EstInfo.cost;
        
        
        sigsqr = EstInfo.sigsqrhat;
        
        
        efit= gof(ytrue,ye);
        gfit = gof(gtrue,ghat);
        Efit = [Efit;efit];
        Cost = [Cost;cost];
        Gfit = [Gfit;gfit];
        
        Sigsqr = [Sigsqr;sigsqr];
        
        Bv = zeros(dN-n+1,M);
        for i = 1:M
            Bv(:,i) = uv(n:dN).^i;
        end
        bv = Bv*poly;
        % yp = Ov*W + hyper(end) + bv;
        yp = Ov*W + bv;
        pfit = gof(yv_true,yp);
        Pfit = [Pfit;pfit];
        Hp = [Hp hyper];
        Poly = [Poly poly];
        rg = (-10:0.1:10)';
        poly_hat = [0;poly]; %poly_hat = [hyper(end);poly];
        ftrue = zeros(length(rg),1);
        fhat = zeros(length(rg),1);
        for i = 0:M
            ftrue = ftrue + poly_true(i+1)*rg.^i;
            fhat = fhat + poly_hat(i+1)*rg.^i;
        end
        plfit = gof(ftrue, fhat);
        Plfit = [Plfit;plfit];
        
        
    end
    fprintf('Completed\n');
    toc;
    
    EFIT = [EFIT Efit];
    PFIT = [PFIT Pfit];
    GFIT = [GFIT Gfit];
    PLFIT = [PLFIT Plfit];
    COST = [COST Cost];
    SIGSQR = [SIGSQR Sigsqr];
end

p_boxplot(EFIT',0,100,kernel,'Efit','');

p_boxplot(PFIT',0,100,kernel,'Pfit','');

p_boxplot(GFIT',0,100,kernel,'Gfit','');

p_boxplot(PLFIT',0,100,kernel,'Plfit','');
save('Results/EFIT.mat','EFIT');
save('Results/PFIT.mat','PFIT');
save('Results/GFIT.mat','GFIT');
save('Results/PLFIT.mat','PLFIT');
save('Results/COST.mat','COST');
save('Results/Hp.mat','Hp');

pp = zeros(M,Maxrepi,length(kernel));
for i = 1:Maxrepi
    for j = 1:length(kernel)
        pp(:,i,j) = Poly(:,(i-1)*length(kernel)+j);
    end
end
PP = [pp(:,:,1)' pp(:,:,2)' pp(:,:,3)'];

save('Results/PP.mat', 'PP');
% p_boxplot(PP,-5,5,{'DC_bd:p1','DC_bd:p2','DC_opt:p1', 'DC_opt:p2', 'DC_dc:p1', 'DC_dc:p2',},'Poly','');
p_boxplot(PP,-5,5,{'DC_bd:p1','DC_bd:p2', 'DC_bd:p3','DC_opt:p1', 'DC_opt:p2','DC_opt:p3', 'DC_dc:p1', 'DC_dc:p2', 'DC_dc:p3'},'Poly','');
%p_boxplot(PP,-5,5,{'DC_bd:p1','DC_bd:p2','DC_opt:p1', 'DC_opt:p2', 'DC_dc:p1', 'DC_dc:p2',...
%     'NEW_bd:p1','NEW_bd:p2','NEW_new:p1','New_new:p2',},'Poly','');
% p_boxplot(PP,-5,5,{'NEW_bd:p1','NEW_bd:p2','NEW_new:p1','New_new:p2',},'Poly','');

save('Results/SIGSQR_gmleC.mat','SIGSQR')