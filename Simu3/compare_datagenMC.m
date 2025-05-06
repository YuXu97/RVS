clear;

addpath('generatedata/');
N= 250;
par.finalTime       = N;               % Data length
par.resampling      = 3;               % 1 - multinom, 2 - systematic, 3 - stratified
par.Np              = 15;              % Number of particles
par.numMCMC         = 20000;           % Number of MCMC iterations
par.burnin = 10000;                     % MCMC burnin. Also, stops rescaling after this many iterations.
par.input           = 1;               % Use this to indicate whether or not we should use an input signal when generating data
par.scaling         = 'h2';            % How to resolve the scale ambiguity, either 'range', 'H2' or 'none'

%
par.dynprior        = 'ard';
par.ardnx           = 15;


Maxrepi = 40;
parfor repi = 1:Maxrepi
    testsystem_genMC(N, par, repi);
end