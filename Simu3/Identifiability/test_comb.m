% clear;
% SIGSQR_gmleC = SIGSQR;
ii = 3;
addpath([pwd '/BasicFunctions']);
addpath([pwd '/Databank']);
sigsqr_true = zeros(80,1);
for i = 1:80
    d = load(['Databank/data_N3000_repi=' int2str(i) '.mat']);
    sigsqr_true(i) = d.datainfo.NoiseVariance;
end
  
plot(sigsqr_true); hold on;
plot(SIGSQR_ml(ii,:)); hold on;
plot(SIGSQR_gmle(ii,:)); hold on;
plot(SIGSQR_gmleC(ii,:)); grid on;
legend('true','ml','gmle','gmleC');

gof(sigsqr_true,SIGSQR_ml(ii,:))
gof(sigsqr_true,SIGSQR_gmle(ii,:))
gof(sigsqr_true,SIGSQR_gmleC(ii,:))