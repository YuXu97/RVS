
%s0ve('Results/0HK_M3/ind.m0t','ind');
% 
% EFITc = EFIT;
% PFITc = PFIT;
% GFITc = GFIT;
% PLFITc = PLFIT;
% COSTc = COST;COST
% HPc = Hp;
% PPc = PP;
% 
% save('Results/CHK_M2/EFITc.mat','EFITc');
% save('Results/CHK_M2/PFITc.mat','PFITc');
% save('Results/CHK_M2/GFITc.mat','GFITc');
% save('Results/CHK_M2/PLFITc.mat','PLFITc');
% save('Results/CHK_M2/COSTc.mat','COSTc');
% save('Results/CHK_M2/HPc.mat','HPc');
% save('Results/CHK_M2/PPc.mat','PPc');
% 
% % 
% plot(1:20,COSTa(3,:)); hold on;
% plot(1:20,COSTc(3,:)); grid on;
% 
% legend('3','500')
figure(1);
plot(COST0(3,:)); hold on;
plot(COSTa(3,:)); hold on;
plot(COSTc(3,:)); grid on;
legend('ratio=3','ratio=10','ratio=500');
title('DC-dc: COST');

figure(2);
subplot(1,2,1);
plot(EFIT0(3,:)); hold on;
plot(EFITa(3,:)); hold on;
plot(EFITc(3,:)); grid on;
legend('ratio=3','ratio=10','ratio=500');
title('DC-dc: Efit');


subplot(1,2,2);
plot(PFIT0(3,:)); hold on;
plot(PFITa(3,:)); hold on;
plot(PFITc(3,:)); grid on;
legend('ratio=3','ratio=10','ratio=500');
title('DC-dc: Pfit');

figure(3);
subplot(1,2,1);
plot(GFIT0(3,:)); hold on;
plot(GFITa(3,:)); hold on;
plot(GFITc(3,:)); grid on;
legend('ratio=3','ratio=10','ratio=500');
title('DC-dc: Gfit');

subplot(1,2,2);
plot(PLFIT0(3,:)); hold on;
plot(PLFITa(3,:)); hold on;
plot(PLFITc(3,:)); grid on;
legend('ratio=3','ratio=10','ratio=500');
title('DC-dc: Polyfit');

