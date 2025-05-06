function [A, B, C, h] = identify_wienergp(sys0,data,par)

model.type = sys0.type;
ard = strcmpi(par.dynprior,'ard');
h2scaling = strcmpi(par.scaling,'h2');
rangescaling = strcmpi(par.scaling,'range');
nu = size(sys0.ss.B,2);

if(par.save)
    savethese = {'data','sys0','par','k','x','J','z','A','Q','h','R','np','zp','lambda','zeta','prior','proposal'};
    if(nu > 0)
        savethese = {savethese{:},'B'};
    end
    
    if(ard)
        savethese = {savethese{:},'delta'};
    end
end


%% PRIOR ===================================================================

if(ard)
    nx = par.ardnx; % Maximum order considered
else
    nx = sys0.nx;   % If MNIW prior, nx is assumed known
end
    
model.nx = nx;
model.ss.C = [eye(sys0.ny) zeros(sys0.ny,nx-sys0.ny)]; % C = [1 0 ... 0] always (single output)
model.ss.X1 = zeros(nx,1); % We need to do something with the initial condition!
model.ss.P1 = eye(nx); 

empcov = cov(data.y'); % Compute the empirical covariance of the data
ny = sys0.ny;

fprintf('Initialising using N4SID... \n');
dat = iddata(data.y',data.u');
initSys = n4sid(dat, nx);
if(any(abs(eig(initSys.A)) > 1))
    fprintf('Unstable system. Projecting poles');
    [V,D] = eig(initSys.A);
    D(abs(D) > 1) = 1./D(abs(D) > 1);
    initSys.A = real(V*D/V);
end
initSys = ss(initSys.A, initSys.B, initSys.C, zeros(1,nu));
initSys = obscanon(initSys);
fprintf('Done!\n');

if(ard)
    % ARD for [A B] (i.e. gamma for delta)
    %prior.a = nx; 
    prior.a = nx*ones(1,nx+nu);
    prior.b = nx/1000;
    prior.a(1) = prior.b/10;        % First state component is always included
    prior.a(nx+1:end) = prior.b/10; % All inputs are always included
else
    % MN for [A B]
    prior.M = [initSys.a initSys.b];
    prior.L = eye(nx+nu);
end
% IW for Q
prior.n0 = nx+2;
S0 = 0.625*empcov; % See Fox (2009) p. 160
prior.S0 = blkdiag(S0, (det(S0)^(1/(nx-ny)))*eye(nx-ny));
% IW for R
prior.r0 = ny+2; % See Fox (2009) p. 160
prior.R0 = 0.075*empcov;
% GP for h(.) - we use the GPML syntax
prior.zetaSigma = 1000000*eye(2); % Covariance for GP hyperprior 
proposal.zetaSigma = 0.1*eye(2); % Covariance for random walk MH kernel for GP hyperparameters
meanfunc = @meanLinear;
hyp.mean = 1; % linear prior mean 
%covfunc = @my_covEiso; % isotropic exponential covariance
covfunc = @covSEiso; % isotropic Gaussian covariance

%% STEP 1: Initialise ======================================================
% Number of points at which we evaluate the GP - these are only used to
% visualise the posterior of h
np = 200;
zp = linspace(3*floor(min(data.z)), 3*ceil(max(data.z)+1), np);
h0 = sys0.ss.h(zp);
lambda = max(data.z)-min(data.z); % Used to scale the system (if range-scaling is enabled)

T = par.finalTime;
% Allocate memory
A = zeros(nx,nx,par.numMCMC);
B = zeros(nx,nu,par.numMCMC);
Q = zeros(nx,nx,par.numMCMC);
if(ard), delta = zeros(1,nx+nu,par.numMCMC); end;
R = zeros(ny,ny,par.numMCMC);
h = zeros(1,np,par.numMCMC);
zeta = zeros(1,2,par.numMCMC);
LYN = zeros(1,par.numMCMC);

% Initialise parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GP hyperparameters
ell = 5; % lengthscale
sf = 5;  % amplitude
zetatmp = log([ell  sf]);

A(:,:,1) = initSys.A;
Btmp = initSys.B;
Qtmp = prior.S0;
R(:,:,1) = prior.R0;
htmp = zp;

model.ss.A = A(:,:,1);
model.ss.B = Btmp;
model.ss.Q = Qtmp;
model.ss.R = R(:,:,1);
model.ss.h = @(z)(gpinterp(z,zp,htmp));

% Initialise state trajectory
[x,J,tmpLYN] = bsi_wiener(data, model, par);
LYN(1) = sum(tmpLYN);

% Resolve scale ambiguity
z = model.ss.C*x;   % Input to nonlinearity for fixed state trajectory
if(h2scaling) % H2 norm
    H2 = cmpH2norm(A(:,:,1), Btmp, model.ss.C, Qtmp, 0);
    sc = sys0.H2/H2;
elseif(rangescaling) % Range
    sc = lambda/(max(z)-min(z)); % 
else
    sc = 1;
end
x = sc*x;
z = sc*z;
B(:,:,1) = sc*Btmp;
Q(:,:,1) = sc^2*Qtmp;
h(:,:,1) = model.ss.h(zp/sc); % Reevaluate the function so that "htilde" is evaluated at zp
zeta(:,:,1) = zetatmp + [log(sc) 0]; % Length scale is multiplied by sc
hyp.cov = zeta(:,:,1);

% PLOT
if(par.plotON)
    [mag0, phase0, w] = initplot(sys0, A(:,:,1), B(:,:,1), model.ss.C, Q(:,:,1), R(:,:,1), zp, h0, h(:,:,1), par);
    figure(4);
    if(ard)
        figure(5);
        lastfig = 5;
    else
        lastfig = 4;
    end
    if(par.input)
        spreadfigures(1:lastfig);
    else
        spreadfigures(2:lastfig);
    end
end

%% STEP 2: PMCMC ===========================================================
CC = 1;
CCC = 1;
H = {};
for(k = 2:par.numMCMC)   
    % R - references to Fox (2009) ----------------------------------------
    % Start with R - this value is used when sampling h later on
    hz = model.ss.h(z/sc);
    e = data.y - hz;
    SR = e*e';                                  % (4.21)
    R(:,:,k) = iwishrnd(SR+prior.R0, T+prior.r0);
    
    if(ard) % ARD
        % DRAW A/B | x_{1:T}, Q, delta ------------------------------------
        psi = x(:,2:end);
        psibar = [x(:,1:end-1) ; data.u(1:end-1)];
        ZbarTQ = kron(psibar, inv(Q(:,:,k-1)));
        Zbar = kron(psibar', eye(nx));
        Sig0i = kron(diag(delta(1,:,k-1)), eye(nx));
        SigA = inv(Sig0i + ZbarTQ*Zbar);
        SigA = (SigA+SigA')/2;
        muA = SigA*ZbarTQ*psi(:);
        
        vecAB = mvnrnd(muA', SigA, 1)';
        AB = reshape(vecAB,nx,nx+nu);
        A(:,:,k) = AB(:,1:nx);
        Btmp = AB(:,nx+1:end);
        
        % DRAW delta | A/B
        a = prior.a + nx/2;
        b = prior.b + 1/2*sum(AB.^2,1);
        delta(1,:,k) = gamrnd(a,1./b);
        
        % DRAW Q | x_{1:T}, A/B
        err = psi - AB*psibar;
        S4 = err*err';
        Qtmp = iwishrnd(S4+prior.S0, T-1+prior.n0); % (4.11)
    else % MNIW
        % DRAW {A/B, Q} | x_{1:T} -----------------------------------------
        psi = x(:,2:end);
        psibar = [x(:,1:end-1) ; data.u(1:end-1)];
        S1 = psibar*psibar' + prior.L;              % (4.10)
        S2 = psi*psibar' + prior.M*prior.L;
        S3 = psi*psi' + prior.M*prior.L*prior.M';
        L = chol(S1); V = S2/L;
        % S4 = S3 - (S2/S1)*S2'
        S4 = S3-V*V';                               % Below (4.11)
        Qtmp = iwishrnd(S4+prior.S0, T-1+prior.n0); % (4.11)
        AB = mxnormrnd(S2/S1, Qtmp, S1);      % (4.9)
        A(:,:,k) = AB(:,1:nx);
        Btmp = AB(:,nx+1:end);
    end
    
    % zeta | r, x_{1:T} ---------------------------------------------------
    zeta_star = zeta(:,:,k-1) + mvnrnd([0 0], proposal.zetaSigma);    
    zval = z';
    m = feval(meanfunc, hyp.mean, zval);
    eval = data.y'-m;
    
    K0 = feval(covfunc, zeta(:,:,k-1), zval) + R(:,:,k)*eye(T);
    K0 = (K0+K0')/2;
    K1 = feval(covfunc, zeta_star, zval) + R(:,:,k)*eye(T);
    K1 = (K1+K1')/2;
    
    R0 = chol(K0);
    R1 = chol(K1);
    
    s0 = R0'\eval;
    s1 = R1'\eval;
    prob = det(R1\R0)*exp(1/2*(s0'*s0)-1/2*(s1'*s1));
    
    accept = (rand(1) < prob);
    if(accept)
        zetatmp = zeta_star;        
        hyp.cov = zetatmp;
    else
        zetatmp = zeta(:,:,k-1);
        hyp.cov = zetatmp; % Set hyp.cov here, since we might have rescaled zeta
    end
    
    
    % h(.) - Gaussian process posterior -----------------------------------
    hyp.lik = log(sqrt(R(:,:,k)));
    % Note that we have rescaled z from the previous iteration
    [mu, sqS] = mygp(hyp, meanfunc, covfunc, z', data.y', zp');
    %htmp = mvnrnd(mu,S)';
    htmp = mu + sqS'*randn(np,1);
    
    % Run a filter and sample a new trajectory ----------------------------
    % Parameterise the model with the updated parameters
    model.ss.A = A(:,:,k);
    model.ss.B = Btmp;
    model.ss.Q = Qtmp;
    model.ss.R = R(:,:,k);
    model.ss.h = @(z)(gpinterp(z,zp,htmp));
    H{k} = model.ss.h;
    
    [x,J,tmpLYN] = bsi_wiener(data, model, par, x, J);
    LYN(k) = sum(tmpLYN);
    
    % Resolve scale ambiguity ---------------------------------------------
    z = model.ss.C*x;   % Input to nonlinearity for fixed state trajectory
    if(k <= par.burnin)
        if(h2scaling) % H2 norm
            H2 = cmpH2norm(A(:,:,k), Btmp, model.ss.C, Qtmp, 0);
            sc = sys0.H2/H2;
        elseif(rangescaling)
            sc = lambda/(max(z)-min(z)); %
        else
            sc = 1;
        end
    else    
        % Only rescale during burnin phase to make sure that the Markov chain
        % gets the correct stationary distribution        
        sc = 1;
    end
        
    x = sc*x;
    z = sc*z;
    B(:,:,k) = sc*Btmp;
    Q(:,:,k) = sc^2*Qtmp;
    h(:,:,k) = model.ss.h(zp/sc); % Reevaluate the function so that "htilde" is evaluated at zp
    zeta(:,:,k) = zetatmp + [log(sc) 0]; % Length scale is multiplied by sc
        
    if(par.plotON),
        updateplot(k, A(:,:,k), B(:,:,k), model.ss.C, Q(:,:,k), R(:,:,k), zp, h0, h(:,:,k), z, data, zeta(:,:,k), mag0, phase0, w, par);
        if(ard) % Plot delta
            figure(5);
            burnin = max(min(5000,k-1000),0);
            boxplot(squeeze(delta(1,:,burnin+1:k))');
            title(sprintf('ARD precisions - using samples %i to %i',burnin+1,k));
            drawnow;
        end
    end;
    
    if(k >= CC*100)
        fprintf('%i: R/R0 = %f\n',k,R(:,:,k)/sys0.ss.R);
        CC = CC + 1;
    end
    if(par.save && k >= CCC*1000)
        save(par.filename,savethese{:});
        CCC = CCC + 1;
    end
    
end

A = A(:,:,par.burnin+1:end);
B = B(:,:,par.burnin+1:end);
% A = mean(A(:,:,par.burnin+1:end),3);
% B = mean(B(:,:,par.burnin+1:end),3);
C = model.ss.C;
h = H(par.burnin+1:end);
% h = mean(h(:,:,par.burnin+1:end),3)
end

%--------------------------------------------------------------------------
function updateplot(k, A, B, C, Q, R, zp, h0, h, z, data, zeta, mag0, phase0, w, par)
% BODE PLOT -------------------------------------------------------
if(par.input)
    figure(1);
    [mag, phase] = bode(ss(A,B,C,0,1),w);
    subplot(2,1,1);
    plot(w(:), 20*log10(mag0(:)),'k-');
    hold on;
    plot(w(:), 20*log10(mag(:)),'b-');
    xlim([0, pi]);
    hold off;
    title('Bode plot / amplitude');
    subplot(2,1,2);
    plot(w(:), phase0(:),'k-');
    hold on;
    plot(w(:), phase(:),'b-');
    xlim([0, pi]);
    hold off;
    title('Phase');
    drawnow
end
% PLOT NOISE COVARIANCE DETERMINANTS ------------------------------
figure(2);
subplot(2,1,1);
%plot(k,det(Q),'k.');
semilogy(k,sqrt(trace(Q)),'k.');
xlim([0,100*ceil(k/100)]);
subplot(2,1,2);
semilogy(k,sqrt(R),'k.');
xlim([0,100*ceil(k/100)]);
title('sqrt(R)');
drawnow;
% PLOT NONLINEAR MAPPING ------------------------------------------
figure(3);
plot(zp,h0,'k-');
hold on;
plot(zp,h,'b-');
plot(z,data.y,'b.');
plot(data.z, data.y,'k.');
xlim([floor(min(z)), ceil(max(z))]);
hold off;
title('Nonlinearity h');

drawnow;
% PLOT HYPERPARAMETERS AND SCALING --------------------------------
figure(4);
subplot(211);
plot(k,exp(zeta(1)),'k.'); hold on;
title('GP kernel length scale');
subplot(212);
plot(k,exp(zeta(2)),'k.'); hold on;
title('GP kernel amplitude');
end

%--------------------------------------------------------------------------
function [mag0, phase0, w] = initplot(sys0, A, B, C, Q, R, zp, h0, h, par)
% BODE PLOT ----------------------------------------------------------
if(par.input)
    figure(1);clf(1);
    [mag0, phase0, w] = bode(ss(sys0.ss.A,sys0.ss.B,sys0.ss.C,0,1));
    [mag, phase] = bode(ss(A,B,C,0,1),w);
    subplot(2,1,1);
    plot(w(:), 20*log10(mag0(:)),'k-');
    hold on;
    plot(w(:), 20*log10(mag(:)),'b-');
    xlim([0, pi]);
    hold off;
    title('Bode plot / amplitude');
    subplot(2,1,2);
    plot(w(:), phase0(:),'k-');
    hold on;
    plot(w(:), phase(:),'b-');
    xlim([0, pi]);
    hold off;
    title('Phase');
    drawnow
end
% PLOT NOISE COVARIANCE DETERMINANTS/TRACES ---------------------------
figure(2);clf(2);
subplot(2,1,1);
semilogy([0, par.numMCMC],sqrt(trace(sys0.ss.Q))*[1 1],'k-');
hold on;
semilogy(1,trace(Q),'k.');
xlim([0,100]);
title('sqrt(trace(Q))');

subplot(2,1,2);
semilogy([0, par.numMCMC],sqrt(sys0.ss.R)*[1 1],'k-');
hold on;
plot(1,R,'k.');
xlim([0,100]);
title('sqrt(R)');

% PLOT NONLINEAR MAPPING ----------------------------------------------
figure(3);
plot(zp,h0,'k-');
hold on;
plot(zp,h,'b-');
hold off;
title('Nonlinearity h');
end

