clear;
Maxrepi = 40;
% c=parcluster('local');
% c.NumWorkers= 40;
% parpool(40);

addpath('generatedata/')
addpath('pmcmc/')
addpath('helpers/')
addpath(genpath('gpml-matlab-v3.1-2010-09-27'))

n = 100;
N = 250;

par.finalTime       = N;            % Data length
par.resampling      = 3;               % 1 - multinom, 2 - systematic, 3 - stratified
par.Np              = 15;              % Number of particles
par.numMCMC         = 20000;           % Number of MCMC iterations
par.burnin = 10000;                     % MCMC burnin. Also, stops rescaling after this many iterations.
par.input           = 1;               % Use this to indicate whether or not we should use an input signal when generating data
par.scaling         = 'h2'; %'h2';            % How to resolve the scale ambiguity, either 'range', 'H2' or 'none'
%
par.dynprior        = 'ard';

par.plotON      = 0;
par.save        = 0;

PFIT_gp = [];
% list = [5;7;16;17;33;35;37];
parfor repi = 1:Maxrepi  
%   kkk = 1:length(list)%
%     repi = list(kkk);

    d = load(['databank/' 'data_' int2str(repi) '.mat']);
    s = load(['databank/' 'sys0_' int2str(repi) '.mat']);
    s.sys0.H2 = cmpH2norm(s.sys0.ss.A, s.sys0.ss.B, s.sys0.ss.C, s.sys0.ss.Q, 0);

%     order = 10;
%     err = zeros(length(order),1);
%     for itero = 1:length(order)
%         par_iter = par;
%         par_iter.ardnx  = order(itero);
%         [A,B,C,h] = identify_wienergp(s.sys0,d.data, par_iter);
% 
%         yhat = zeros(2*N,1);
%         for iter = 1:par.numMCMC - par.burnin
%             Ai = A(:,:,iter);
%             Bi = B(:,:,iter);
%             ht = h(iter);
%             htt = ht{:};
% 
%             sys_est = ss(Ai,Bi,C,0,1);
% 
%             z = sim(idss(sys_est), [d.data.u'; d.data.uv']);
% 
%             warning off
%             yhat = yhat + htt(z)';
%         end
%         yhat = yhat/(par.numMCMC - par.burnin);
%         err(iter) = sumsqr(yhat(N+1:2*N)-d.data.yv');
%     end
%     gof(d.data.yvtrue(n:end)',yhat(N+n:2*N))
% % 
%     [~, min_id] = min(err(:));
%     order_opt = order(min_id);

    order_opt = 10;
    par_opt = par;
    par_opt.finalTime = 2*N;
    par_opt.ardnx  = order_opt;
    d_opt = d;
    d_opt.data.u = [d.data.u, d.data.uv];
    d_opt.data.y = [d.data.y, d.data.yv];
    [A,B,C,h] = identify_wienergp(s.sys0, d_opt.data, par_opt);

    yhat = zeros(4*N,1);
    ghat = zeros(n,1);
    xx = (-1.5:1e-2:1.5)';
    phat = zeros(length(xx),1);

    gtrue = d.data.gtrue(1:n);
    for iter = 1:par.numMCMC - par.burnin
        Ai = A(:,:,iter);
        Bi = B(:,:,iter);
        ht = h(iter);
        htt = ht{:};

        sys_est = ss(Ai,Bi,C,0,1);
        
        g = impulse(sys_est,0:n-1);

        sc = gtrue(2)/g(2);
        
        z = lsim(sys_est, [d_opt.data.u'; d_opt.data.ut'], 1:4*N);
        p = htt(xx/sc)';

        
        warning off
        ghat = ghat + g*sc;
        yhat = yhat + htt(z)';
        phat = phat + p;
    end
    yhat = yhat/(par.numMCMC - par.burnin);
    ghat = ghat/(par.numMCMC - par.burnin);
    phat = phat/(par.numMCMC - par.burnin);
    
    % PFIT_gp = [PFIT_gp gof(d.data.yttrue(n:end)',yhat(2*N+n:4*N))];
    
    
    gfit = gof(gtrue,ghat);
    nfit = gof(d.data.nonl(xx),phat);
    pfit = gof(d.data.yttrue(n:end)',yhat(2*N+n:4*N));
    s = struct('pfit',pfit);
    save(sprintf('Results/6dsys/gp/PFIT_gp_%d.mat',repi),'-fromstruct',s);
    sg = struct('gfit',gfit);
    save(sprintf('Results/6dsys/gp/GFIT_gp_%d.mat',repi),'-fromstruct',sg);
    sn = struct('nfit',nfit);
    save(sprintf('Results/6dsys/gp/NFIT_gp_%d.mat',repi),'-fromstruct',sn);
end
% save('Results/4dsys/gp/PFIT_gp.mat','PFIT_gp');
% save('Results/6dsys/gp/PFIT_gp.mat','PFIT_gp');
