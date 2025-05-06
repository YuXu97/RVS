clear;

addpath('generatedata/')
% c=parcluster('local'); c.NumWorkers= 3; parpool(3);


Maxrepi = 1;
N = 500;
n = 80;

PFIT_nlwh = [];

for repi = 1:Maxrepi
    fprintf('repi = %i: \n',repi);
    d = load(['databank/' 'data_' int2str(repi) '.mat']);
    data = [d.data.u(1:N)',d.data.y(1:N)'];
    data = flip(data,2);

%     %hold out cross validation
    order = 1:40;
    err = zeros(length(order),1);
    for iter = 1:length(order)
        ord = order(iter);
        ORDs = [ord ord 0];
        sys = nlhw(data, ORDs, [], idSaturation);
        y_nlhw = sim(sys, [d.data.u'; d.data.uv']);
        err(iter) = sumsqr(y_nlhw(N+1:2*N) - d.data.yv');
        %gof(d.data.yvtrue, y_nlhw)
    end

    [~, min_id] = min(err(:));
    order_opt = order(min_id);

    ORD_opt = [order_opt,order_opt,1];
    data_est = [[d.data.u d.data.uv]',[d.data.y d.data.yv]'];
    data_est = flip(data_est,2);
    sys_opt = nlhw(data_est, ORD_opt, [], idSaturation);
    y_nlhw_opt = sim(sys_opt, [d.data.u d.data.uv d.data.ut]');
    PFIT_nlwh = [PFIT_nlwh, gof(d.data.yttrue(n:end), y_nlhw_opt(2*N+n:4*N))];
    fprintf('Completed\n');
end

save('Results/nlwh/PFIT_nlwh.mat','PFIT_nlwh');