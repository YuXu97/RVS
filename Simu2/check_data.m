clear;
Nmax = 4000;
Maxrepi = 200;
lack_list = [];
for i = 1:Maxrepi
    try 
        d = load(['Databank/data_N' int2str(Nmax) '_repi=' int2str(i) '.mat']);
    catch
        lack_list = [lack_list;i];
    end
    
end