
lb = 0;
ub = 100;

Xe = [Efit_DC_bd Efit_DC_opt Efit_DC_opt_dc];
Xp = [Pfit_DC_bd Pfit_DC_opt Pfit_DC_opt_dc];

[~,columns] = size(Xe);

[Length,~] = size(Efit_DC_bd);
names = [ repmat({'DC-bd'},Length,1); repmat({'DC-opt'},Length,1); repmat({'DC-dc'},Length,1)];

names = categorical(names);

% boxchart(names,Xe)

% 
figure('position', [500, 500, 800, 600]) 
tiledlayout(1,2)



% Left axes
ax1 = nexttile;
b1 = boxchart(ax1,Xe(:),'GroupByColor',names,'MarkerStyle','.', 'BoxWidth', 0.7);
ylim([lb ub]);
ylabel(ax1,'Estimation Fit')
legend
grid on;
xlabel('Training (80)')
ylim1=get(gca,'ylim');
% set(gca,'YTick',93:0.5:100)
n_outliers = sum(Xe < lb);
for i = 1 : columns
    temp_outliers = n_outliers(i);
    if temp_outliers ~= 0
    str_temp = ['+', int2str(temp_outliers)];
    yvalue = lb + diff(ylim1)*.03;
    text(0.4+0.32*i,yvalue,str_temp);
    end
end


% Right axes
ax2 = nexttile;
b2 = boxchart(ax2,Xp(:),'GroupByColor',names,'MarkerStyle','.', 'BoxWidth', 0.7);
ylim([lb ub]);
ylabel(ax2,'Predictiontion Fit')
legend
grid on;
xlabel('Validation (560)')
ylim2=get(gca,'ylim');
% set(gca,'YTick',93:0.5:100)
n_outliers = sum(Xp < lb);
for i = 1 : columns
    temp_outliers = n_outliers(i);
    if temp_outliers ~= 0
    str_temp = ['+', int2str(temp_outliers)];
    yvalue = lb + diff(ylim2)*.03;
    text(0.4+0.32*i,yvalue,str_temp);
    end
end


