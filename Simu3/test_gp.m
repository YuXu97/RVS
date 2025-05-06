clear;
%
% c=parcluster('local');
% c.NumWorkers= 75;
% parpool(75);


addpath('generatedata/')
addpath('pmcmc/')
addpath('helpers/')
addpath(genpath('gpml-matlab-v3.1-2010-09-27'))
addpath(genpath('No_Identifiability_gmle'))
n = 70;
N = 400;


par.finalTime       = N;            % Data length
par.resampling      = 3;               % 1 - multinom, 2 - systematic, 3 - stratified
par.Np              = 15;              % Number of particles
par.numMCMC         = 20000;           % Number of MCMC iterations
par.burnin = 10000;                     % MCMC burnin. Also, stops rescaling after this many iterations.
par.input           = 1;               % Use this to indicate whether or not we should use an input signal when generating data
par.scaling         = 'h2'; %'h2';            % How to resolve the scale ambiguity, either 'range', 'H2' or 'none'

%
par.dynprior        = 'ard';
par.ardnx           = 10;


par.plotON      = 0;
par.save        = 0;

% for repi = 1:1
%     testsystem_wienergp_2d(N, par, repi);
% end
EFIT = [];
PFIT = [];
np = 200;

parfor repi = 1:160
    d = load(['databank/' 'data_' int2str(repi) '.mat']);
    s = load(['databank/' 'sys0_' int2str(repi) '.mat']);
    zp = linspace(3*floor(min(d.data.z(1:N))), 3*ceil(max(d.data.z(1:N))+1), np);

    s.sys0.H2 = cmpH2norm(s.sys0.ss.A, s.sys0.ss.B, s.sys0.ss.C, s.sys0.ss.Q, 0);
    break_flag = 0;
    iter_count = 0;
    while ~break_flag
        try
            [A,B,C,h] = identify_wienergp(s.sys0,d.data,par);

            ytrue_e = d.data.y'-d.data.e(1:N)';
            ytrue_v = d.data.yv'-d.data.e(N+1:2*N)';

            yhat = zeros(2*N,1);
            for iter = 1:par.numMCMC - par.burnin
                Ai = A(:,:,iter);
                Bi = B(:,:,iter);
                ht = h(iter);
                htt = ht{:};

                sys_est = ss(Ai,Bi,C,0,1);
                z = lsim(sys_est, [d.data.u'; d.data.uv'], 1:2*N);
                warning off
                yhat = yhat + htt(z)';
            end
            yhat = yhat/(par.numMCMC - par.burnin);

            %         yhat = h(z);
            efit = gof(ytrue_e(n:N),yhat(n:N))
            pfit = gof(ytrue_v(n:N),yhat(N+n:2*N))

            EFIT = [EFIT; efit];
            PFIT = [PFIT; pfit];
            break_flag = 1;
        catch
            break_flag = 0;
        end

        iter_count = iter_count + 1;
        if iter_count > 50 && ~break_flag
            efit = nan;
            pfit = nan;
            EFIT = [EFIT; efit];
            PFIT = [PFIT; pfit];
            break;
        end
    end
end
save('EFIT_gp.mat','EFIT');
save('PFIT_gp.mat','PFIT');

% Al = A(:,:,end);
% Bl = B(:,:,end);
% hl = h(end);
% hll = hl{:};
% sys_est = ss(Al,Bl,C,0,1);
% z = lsim(sys_est, [d.data.u'; d.data.uv'], 1:2*N);
% yhat = hll(z)';
% efit = gof(ytrue_e(n:N),yhat(n:N))
% pfit = gof(ytrue_v(n:N),yhat(N+n:2*N))


%
%     figure(1)
%     plot(ytrue_e(n:N)); hold on;
%     plot(yhat(n:N));grid on;
%     legend('TRUE','ESTIMATE')
%     title('efit')
%
%     figure(2)
%     plot(ytrue_v(n:N)); hold on;
%     plot(yhat(N+n:2*N));grid on;
%     legend('TRUE','ESTIMATE')
%     title('pfit')
% p_boxplot(EFIT(1:80),0,100,'GP','EFIT','')
% p_boxplot(PFIT(1:80),0,100,'GP','PFIT','')