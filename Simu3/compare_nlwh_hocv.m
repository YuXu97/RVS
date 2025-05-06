clear;

addpath('generatedata/')
% c=parcluster('local');
% c.NumWorkers= 3;
% parpool(3);


Maxrepi = 40;
N = 250;
n = 100;

PFIT_nlwh = [];
GFIT_nlwh = [];
NFIT_nlwh = [];
parfor repi = 1:Maxrepi
    fprintf('repi = %i: \n',repi);
    d = load(['databank/' 'data_' int2str(repi) '.mat']);
    data = [d.data.u(1:N)',d.data.y(1:N)'];
    data = flip(data,2);
    gtrue = d.data.gtrue(1:n);

    %hold out cross validation
    order = 1:10;
    order_pl = 9;
    err = zeros(length(order), length(order_pl));
    for iter = 1:length(order)
        ord = order(iter);
        for iter_pl = 1:length(order_pl)
            ord_pl = order_pl(iter_pl);

            ORDs = [ord ord 1];
            pl = idPolynomial1D('Degree', ord_pl);
            sys = nlhw(data, ORDs, [], pl);
            y_nlhw = sim(sys, [d.data.u'; d.data.uv']);
            err(iter,iter_pl) = sumsqr(y_nlhw(N+1:2*N) - d.data.yv');
            %gof(d.data.yvtrue, y_nlhw)
        end 
    end

    [~, min_id] = min(err(:));
    [id_opt, id_opt_pl] = ind2sub([iter,iter_pl], min_id);
    order_opt = order(id_opt);
    order_opt_pl = order_pl(id_opt_pl);

    pl_opt =idPolynomial1D('Degree',order_opt_pl);
    ORD_opt = [order_opt,order_opt,1];

    data_est = [[d.data.u d.data.uv]',[d.data.y d.data.yv]'];
    data_est = flip(data_est,2);
    sys_opt = nlhw(data_est, ORD_opt, [], pl_opt);

    Lsys = sys_opt.LinearModel;
    ghat = impulse(Lsys,0:n-1);
    sc = gtrue(2)/ghat(2);
    ghat = ghat*sc;
    xx = (-1.5:1e-2:1.5)';
    phat = polyval(sys_opt.OutputNonlinearity.Coefficients,xx/sc);
    
    y_nlhw_opt = sim(sys_opt, [d.data.u d.data.uv d.data.ut]');

    gfit = gof(gtrue,ghat);
    nfit = gof(d.data.nonl(xx),phat);

    GFIT_nlwh = [GFIT_nlwh; gfit];
    NFIT_nlwh = [NFIT_nlwh; nfit];
    PFIT_nlwh = [PFIT_nlwh; gof(d.data.yttrue(n:end), y_nlhw_opt(2*N+n:4*N))];
    fprintf('Completed\n');
end

save('Results/6dsys/nlwh/PFIT_nlwh.mat','PFIT_nlwh');
save('Results/6dsys/nlwh/GFIT_nlwh.mat','GFIT_nlwh');
save('Results/6dsys/nlwh/NFIT_nlwh.mat','NFIT_nlwh');