
Maxrepi = 50;
kns = {'DC-bd-w','DC-decay-w','DC-ob-w'};
PFIT = zeros(Maxrepi,length(kns));

for repi = 1:Maxrepi
    Pfit = load([pwd '/Results/Pfit/Pfit' int2str(repi) '.mat']).Pfit;
    pfit = Pfit(:,1)';
    PFIT(repi,:) = pfit;
end

p_boxplot(PFIT,0,100,kns,'','');

