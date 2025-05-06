function [EstInfo] = rvs_ml_smsp(Pi, Rho, y,  n, M, kernel)
[lb,ub] = lubounds(kernel, M);
y = y(n:end);

p = 1;
[~,r] = size(Pi);
pb = p + r;

indm = zeros(M,1);
for m = 1:M
    indm(m) = factorial(pb+m-1)/factorial(m)/factorial(pb-1);
end
indcum = cumsum(indm); indcum = [0;indcum];
CM = zeros(sum(indm),pb);
Cf = zeros(sum(indm),1);
for m = 1:M
    C = MultiPermuations(m,pb);
    CM(indcum(m)+1:indcum(m+1),:) = C;
    for i = 1:indm(m)
        Cf(indcum(m)+i) = multinomial(m,C(i,:));
    end
end


tic;
hpini = ini_rvs(Pi, Rho, y, kernel, M, n, CM, Cf, indcum);
time_ini = toc;
ff = @(x)nglglklhd(x, Pi, Rho, y, kernel, M, n, CM, Cf, indcum);
options = optimoptions('fmincon', 'Display', 'off');

if strcmp(kernel,'DC-bd') || strcmp(kernel,'DC-dc')
    hp = fmincon(ff,hpini,[],[],[],[],lb,ub,@nonlcon_dc,options);
else

%     hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);
    options = optimset('Display', 'off');
    hp = fminsearch(ff,hpini,options);
end




[cost, Oiinvy, O, h0, sigsqr] = nglglklhd(hp, Pi, Rho, y, kernel, M, n, CM, Cf, indcum);


yhat = O*Oiinvy+ h0;
EstInfo.Oiinvy = Oiinvy;
EstInfo.hpini = hpini;
EstInfo.hp = hp;
EstInfo.h0 = h0;
EstInfo.sigsqr = sigsqr;
EstInfo.kernel = kernel;
EstInfo.yhat = yhat;
EstInfo.O = O;
EstInfo.cost = cost;
EstInfo.time_ini = time_ini;
end


function [hpini] = ini_rvs(Pi, Rho, Y, kernel, M, n, CM, Cf, indcum)
%cdiv = 5; lamdiv = 5; %for efficiency test
cdiv = 10; lamdiv = 10; %for accuracy test
if strcmp(kernel, 'TC-bd') || strcmp(kernel, 'DC-bd-sp')

    % hp = [c1 ... CM lam]

    c = linspace(-10,3,cdiv); Lc = length(c);
    lam = linspace(0.01,0.99,lamdiv); Llam = length(lam);
    obj = zeros(Lc,Lc,Lc,Llam);
    for nc1 = 1:Lc
        for nc2 = 1:Lc
            for nc3 = 1:Lc
                for nlam = 1:Llam
                    hpstart = [c(nc1);c(nc2);c(nc3);lam(nlam)];
                    obj(nc1,nc2,nc3,nlam) = nglglklhd(hpstart,Pi, Rho,  Y, kernel, M, n, CM, Cf, indcum);
                end
            end
        end
    end
    [~,indx] = min(obj(:));
    [indc1,indc2,indc3,indlam] = ind2sub([nc1 nc2 nc3 nlam],indx);
    hpini = [c(indc1);c(indc2);c(indc3);lam(indlam)];
%     lb = [-10*ones(M,1); 0.1];
%     ub = [5*ones(M,1); 0.9];
%     div = [10*ones(M,1);10];
%     fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
%     [hpini,~,~] = GridSearch(length(lb), lb, ub, lb+(ub-lb).*0.333, div, 1e-3*ones(length(lb),1), 1e-5, 2000, fun);

end


if strcmp(kernel, 'TC-tc') || strcmp(kernel, 'DC-dc-sp')  || strcmp(kernel, 'DC-ob-sp')

    % hp = [c1 ... CM lam]
    cp = linspace(0.001,8,cdiv); Lcp = length(cp);
    c = linspace(-8,8,cdiv); Lc = length(c);
    lam = linspace(0.01,0.99,lamdiv); Llam = length(lam);
    obj = zeros(Lcp,Lc,Lc,Llam);
    for nc1 = 1:Lcp
        for nc2 = 1:Lc
            for nc3 = 1:Lc
                for nlam = 1:Llam
                    hpstart = [cp(nc1);c(nc2);c(nc3);lam(nlam)];
                    obj(nc1,nc2,nc3,nlam) = nglglklhd(hpstart,Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
                end
            end
        end
    end
    [~,indx] = min(obj(:));
    [indc1,indc2,indc3,indlam] = ind2sub([nc1 nc2 nc3 nlam],indx);
    hpini = [cp(indc1);c(indc2);c(indc3);lam(indlam)];
    %      lb = [0.1; -5*ones(M-1,1); 0.1];
    %      ub = [5; 5*ones(M-1,1); 0.9];
    %      div = [10*ones(M,1);10];
    %      fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
    %      [hpini,~,~] = GridSearch(length(lb), lb, ub, lb+(ub-lb).*0.333, div, 1e-3*ones(length(lb),1), 1e-5, 2000, fun);
end

if strcmp(kernel, 'DC-bd') 

    % hp = [c1 ... CM lam rho]
    ratio = 10;
    lb = [-10*ones(M,1);  0.1; 0.1];
    ub = [5*ones(M,1);  0.9; 1];
    fun = @(x)nglglklhd(x, Pi, Rho,  Y, kernel, M, n, CM, Cf, indcum);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',[0.1*ones(M,1);0.5;0.99], 'objective',fun,'lb',lb,'ub',ub,'nonlcon',@nonlcon_dc);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'DC-dc')  || strcmp(kernel, 'DC-ob')

    % hp = [c1 ... CM lam rho]

    ratio = 10;
    lb = [0.1; -5*ones(M-1,1);  0.1; 0.1];
    ub = [5; 5*ones(M-1,1);  0.9; 1];
    fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0', [0.1*ones(M,1);0.5;0.99], 'objective',fun,'lb',lb,'ub',ub,'nonlcon',@nonlcon_dc);
    hpini = run(ms,problem,ratio*length(ub));
end

end

function [obj, Oiinvy, O, h0, sigsqr] = nglglklhd(hyper, Pi, Rho, y, kernel, M, n, CM, Cf, indcum)
warning('off', 'MATLAB:nearlySingularMatrix');

[N,~] = size(Pi); Nd = N-n+1;

[Pp, Qp] = CalculateOutputKernelGenerators(Pi, Rho, M, n, kernel, hyper, CM, Cf, indcum);

[Q,R] = qr(Pp,0);
bUs = Q;
bVs = Qp*R';
s = sqrt(norm(bUs,'fro')/norm(bVs,'fro'));
if isinf(s)
    obj = 1/eps;
    if nargout > 1
        Oiinvy = zeros(Nd,Nd);
        O = zeros(Nd,Nd);
        h0 = 0;
        sigsqr = 0;
    end
    return
else
    bUs = bUs/s;
    bVs = bVs*s;
end
try
    [Wt,c] = egrss_potrf(bUs',bVs',1);
catch
    obj = 1/eps;
    if nargout > 1
        Oiinvy = zeros(Nd,Nd);
        O = zeros(Nd,Nd);
        h0 = 0;
        sigsqr = 0;
    end
    return
end
logdet = 2*sum(reallog(c));

l1 = egrss_trsv(bUs',Wt,c,y);
le = egrss_trsv(bUs',Wt,c,ones(Nd,1));
h0 = l1'*le/sumsqr(le);

y = y - h0;

Liy = egrss_trsv(bUs',Wt,c,y);

obj = Nd - Nd*log(Nd) + Nd*log(sumsqr(Liy)) + logdet;
if nargout > 1
    Oiinvy = egrss_trsv(bUs',Wt,c,Liy,'T');
    O = tril(Pp*Qp')+triu(Qp*Pp',1);
    sigsqr = sumsqr(Liy)/Nd;
end
end


function [c,ceq] = nonlcon_dc(x)
c = x(end-1) - abs(x(end));
ceq = [];
end
