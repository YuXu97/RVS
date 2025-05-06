
PFIT_gp = zeros(40,1);
for repi = 1:40
d = load(['/home/xuyu/mfiles/VolterraSeriesIDv2/pgwienertar/pgwiener/Results/MC_M2/SMP/ARD10/' 'PFIT_gp_' int2str(repi) '.mat']);
PFIT_gp(repi) = d.pfit;
end
save('/home/xuyu/mfiles/VolterraSeriesIDv2/pgwienertar/pgwiener/Results/MC_M2/SMP/ARD10/PFIT_gp.mat','PFIT_gp')
boxplot(PFIT_gp); ylim([0,100])