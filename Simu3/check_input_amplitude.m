
Amax = 0;
Amin = 0;
for i = 1:40
d = load(['databank/final_nodelay_forvs/data_' int2str(i) '.mat']);
u = [d.data.u d.data.uv d.data.ut]';
Amax = max(Amax, max(u))
Amin = min(Amin, min(u))a
end
