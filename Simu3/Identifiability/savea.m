% 
save('Results/CHK_M2/ind.mat','ind');
% 
EFITa = EFIT(:,ind);
PFITa = PFIT(:,ind);
GFITa = GFIT(:,ind);
PLFITa = PLFIT(:,ind);
COSTa = COST(:,ind);
% 
indt = zeros(1,3*20);
i = 1;
for jj = ind
    indt(i) = 3*(jj-1)+1;
    indt(i+1) = 3*(jj-1)+2;
    indt(i+2) = 3*jj;
    i = i +3 ;
end

HPa = Hp(:,indt);
PPa = PP(ind,:);
% 
save('Results/CHK_M2/EFITa.mat','EFITa');
save('Results/CHK_M2/PFITa.mat','PFITa');
save('Results/CHK_M2/GFITa.mat','GFITa');
save('Results/CHK_M2/PLFITa.mat','PLFITa');
save('Results/CHK_M2/COSTa.mat','COSTa');
save('Results/CHK_M2/HPa.mat','HPa');
save('Results/CHK_M2/PPa.mat','PPa');
% 
% % 
% plot(1:20,EFITa(1,:)); hold on;
% plot(1:20,EFITt(1,:)); grid on;
% legend('after','before')

% EFITn = EFIT(:,neg_id);
% PFITn = PFIT(:,neg_id);
% GFITn = GFIT(:,neg_id);
% PLFITn = PLFIT(:,neg_id);
% PPn = PP(neg_id,:);
% COSTn = COST(:,neg_id);
% 
% p_boxplot(EFITn',0,100,kernel,'Efit-neg','');
% 
% p_boxplot(PFITn',0,100,kernel,'Pfit-neg','');
% 
% p_boxplot(GFITn',0,100,kernel,'Gfit-neg','');
% 
% p_boxplot(PLFITn',0,100,kernel,'Plfit-neg','');


%p_boxplot(PP,-5,5,{'DC_bd:p1','DC_bd:p2', 'DC_bd:p3','DC_opt:p1', 'DC_opt:p2','DC_opt:p3', 'DC_dc:p1', 'DC_dc:p2', 'DC_dc:p3'},'Poly','');
% p_boxplot(PPn,-5,5,{'DC_bd:p1','DC_bd:p2','DC_opt:p1', 'DC_opt:p2', 'DC_dc:p1', 'DC_dc:p2',...
%     'NEW_bd:p1','NEW_bd:p2','NEW_new:p1','New_new:p2',},'Poly-neg','');
% 
% save('Results/EFITn.mat','EFITn');
% save('Results/PFITn.mat','PFITn');
% save('Results/GFITn.mat','GFITn');
% save('Results/PLFITn.mat','PLFITn');
% save('Results/COSTn.mat','COSTn');
% save('Results/PPn.mat','PPn');
% hp1 = zeros(160,1);
% hp2 = zeros(160,1);
% hp3 = zeros(160,1);
% for i = 1:160
%     hp1(i) = Hp(3,3*(i-1)+1);
%     hp2(i) = Hp(3,3*(i-1)+2);
%     hp3(i) = Hp(3,3*i);
% end
% boxplot([hp1 hp2 hp3])
