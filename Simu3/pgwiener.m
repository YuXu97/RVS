function pgwiener(experiment,filename)
% PGWIENER Wiener system identification using PG-BS.
%   This is a procedure for identifying Wiener systems using the approach
%   described in
%       
%       F. Lindsten, T. B. Schön and M. I. Jordan, "Data driven Wiener
%       system identification". Automatica (submitted), 2012.
%
%       F. Lindsten, T. B. Schön and M. I. Jordan, "A semiparametric
%       Bayesian approach to Wiener system identification". In Proceedings
%       of the 16th IFAC Symposium on System Identification, Brussels,
%       Belgium, 2012.
%
%   The procedure contains two example systems.
%
%   >> PGWIENER(k)
%
%   for k = 1, 2 runs the procedure for system k and plots the results. To
%   run the code without plotting, and instead save the identification
%   results to file, run
%
%   >> PGWIENER(k, filename)
%
%   Fredrik Lindsten
%   Linköping, 2012-07-09
%   lindsten@isy.liu.se


addpath(genpath('.'));
close all
%====================
%===   SETTINGS   ===
%====================
par.finalTime       = 1000;            % Data length
par.resampling      = 3;               % 1 - multinom, 2 - systematic, 3 - stratified
par.Np              = 15;              % Number of particles
par.numMCMC         = 25000;           % Number of MCMC iterations
par.burnin = 5000;                     % MCMC burnin. Also, stops rescaling after this many iterations.

par.input           = 1;               % Use this to indicate whether or not we should use an input signal when generating data
par.scaling         = 'h2';            % How to resolve the scale ambiguity, either 'range', 'H2' or 'none'

% Save to file or plot?
if(nargin == 2)
    par.plotON      = 0;
    par.save        = 1;
else
    par.plotON      = 1;
    par.save        = 0;
end

% Choose one of the two available examples
if(experiment == 1)
    par.dynprior        = 'mniw';          % Prior on dynamics, either 'mniw' or 'ard'
    [sys0, data] = testsystem_wienergp_5d(par.finalTime, par);
elseif(experiment == 2)
    par.dynprior        = 'ard';          % Prior on dynamics, either 'mniw' or 'ard'
    par.ardnx           = 10;               % Maximum number of states used if ARD prior
    [sys0, data] = testsystem_wienergp_4d(par.finalTime, par);
else
    error('Unknown experiment. Please select 1 or 2.')
end

sys0.H2 = cmpH2norm(sys0.ss.A, sys0.ss.B, sys0.ss.C, sys0.ss.Q, 0);

%========================================
%===   RUN IDENTIFICATION PROCEDURE   ===
%========================================
fprintf('Starting PMCMC identification of Wiener system.\n\n');

if(par.save)
    par.filename = filename;
    fprintf('\n\nWill save data to file %s!\n\n',par.filename);
else
    dbstop if error;
end

identify_wienergp(sys0,data,par);