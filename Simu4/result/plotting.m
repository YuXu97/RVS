% plot(Nrange,[Time_ini_sp;Time_ini_smsp])
Nrange = 1000:1000:8000;
Pfit_bd = [Pfit_Nfull(1,:);  Pfit_smsp(1,:); Pfit_n(1,:); Pfit_sp(1,:)]';
Pfit_decay = [Pfit_Nfull(2,:);  Pfit_smsp(2,:); Pfit_n(2,:); Pfit_sp(2,:)]';
Pfit_ob = [Pfit_Nfull(3,:); Pfit_smsp(3,:); Pfit_n(3,:); Pfit_sp(3,:)]';

subplot(3,1,1);
nice_plot(Nrange,Time_smsp_ini,1);
kns = {'DC-bd-w','DC-decay-w','DC-ob-w'};
legend(kns)
subplot(3,1,2);
nice_plot(Nrange,Time_sp_ini,1);
subplot(3,1,3);
nice_groupboxplot([Pfit_bd Pfit_decay Pfit_ob],80,100,{'CHOL-N','SMSP','CHOL-n','SP'},{'DC-bd-w','DC-decay-w','DC-ob-w'}, 2, '','',1)
