function [EstInfo] = rvs_ml_sp(Pi, Rho, y,  n, M, kernel)
% [lb,ub] = lubounds(kernel, M);
y = y(n:end);
Pi = Pi(n:end,:); Rho = Rho(1:n,:);

[~,r] = size(Pi);
pb = r;
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

% options = optimoptions('fmincon', 'Display','off');
% hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);


options = optimset('Display', 'off');%, 'MaxFunEvals', 3000, 'MaxIter', 3000, 'TolX', 1e-6,'TolFun',1e-6);
hp = fminsearch(ff,hpini,options);


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
  cdiv = 5; lamdiv = 5;
%  ratio = 1500;
if strcmp(kernel, 'TC-bd') || strcmp(kernel, 'DC-bd-sp')

    % hp = [c1 ... CM lam]

    %     lb = [-10*ones(M,1); 0.1];
    %     ub = [5*ones(M,1); 0.9];
    %     div = [10*ones(M,1);10];
    %     fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
    %     [hpini,~,~] = GridSearch(length(lb), lb, ub, lb+(ub-lb).*0.333, div, 1e-3*ones(length(lb),1), 1e-5, 2000, fun);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    c = linspace(-10,3,cdiv); Lc = length(c);
    lam = linspace(0.01,0.99,lamdiv); Llam = length(lam);
    obj = zeros(Lc,Lc,Lc,Llam);
    for nc1 = 1:Lc
        for nc2 = 1:Lc
            for nc3 = 1:Lc
                for nlam = 1:Llam
                    hpstart = [c(nc1);c(nc2);c(nc3);lam(nlam)];
                    obj(nc1,nc2,nc3,nlam) = nglglklhd(hpstart,Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
                end
            end
        end
    end
    [~,indx] = min(obj(:));
    [indc1,indc2,indc3,indlam] = ind2sub([nc1 nc2 nc3 nlam],indx);
    hpini = [c(indc1);c(indc2);c(indc3);lam(indlam)];
    %%%%%%%%%%%%%%%%%%%%%%%%%

%     lb = [-10*ones(M,1);  0.01];
%     ub = [5*ones(M,1);  0.99];
%     fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
%     ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
%     ms.StartPointsToRun='bounds';
%     ms.Display='off';
%     problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
%     hpini = run(ms,problem,ratio*length(ub));

end


if strcmp(kernel, 'TC-tc') || strcmp(kernel, 'DC-dc-sp') || strcmp(kernel, 'DC-ob-sp')

    %hp = [c1 ... CM lam]

    %     lb = [0.1; -5*ones(M-1,1); 0.1];
    %     ub = [5; 5*ones(M-1,1); 0.9];
    %     div = [10*ones(M,1);10];
    %     fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
    %     [hpini,~,~] = GridSearch(length(lb), lb, ub, lb+(ub-lb).*0.333, div, 1e-3*ones(length(lb),1), 1e-5, 2000, fun);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%
%     lb = [0.1; -10*ones(M-1,1); 0.01];
%     ub = [10; 10*ones(M-1,1); 0.99];
%     fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
%     ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
%     ms.StartPointsToRun='bounds';
%     ms.Display='off';
%     problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
%     hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'DC-bd')

    % hp = [c1 ... CM lam rho]
    ratio = 10;
    lb = [-10*ones(M,1);  0.1; 0.1];
    ub = [5*ones(M,1);  0.9; 1];
    fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'DC-dc') || strcmp(kernel, 'DC-ob')

    % hp = [c1 ... CM lam rho]

    ratio = 10;
    lb = [0.1; -5*ones(M-1,1);  0.1; 0.1];
    ub = [5; 5*ones(M-1,1);  0.9; 1];
    fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0', lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end
end

function [obj, Oiinvy, O, h0, sigsqr] = nglglklhd(hyper, Pi, Rho, y, kernel, M, n, CM, Cf, indcum)
warning('off', 'MATLAB:nearlySingularMatrix');

[Nd,~] = size(Pi);

[Pp, Qp] = CalculateOutputKernelGenerators(Pi, Rho, M, n, kernel, hyper, CM, Cf, indcum);

[~,gamma] = size(Pp);
Ig = eye(gamma);
M = Ig+Qp'*Pp;
% if cond(M) <1e4
%     Mi = Ig/M;
% else
%     try
%         L = chol(M)';
%         Linv = Ig/L;
%         Mi = Linv'*Linv;
%     catch
%         try
%             L = chol(1e-4*Ig+M)';
%             Linv = Ig/L;
%             Mi = Linv'*Linv;
%         catch
%             [Ui,Si,Vi] = svd(M);
%             Mi = Vi*diag(1./diag(Si))*Ui';
%         end
%     end
% end


yp = Pp'*y;
pe = Pp'*ones(Nd,1); qe = Qp'*ones(Nd,1);
Miqe = M\qe; 

h0 = (sum(y)-yp'*Miqe)/(Nd-pe'*Miqe);
Yh = y - h0;
ysp = Pp'*Yh; ysq = Qp'*Yh;
Miysq = M\ysq; 

qf = sumsqr(Yh) - ysp'*Miysq;

if qf<0
    obj = 1/eps;
    h0 = 0;
    Oiinvy = zeros(Nd,Nd);
    O = zeros(Nd,Nd);
    sigsqr = 0;
    return
end

ld = logdet(M);
obj = Nd*(1-log(Nd)) + Nd*log(qf) + ld;

if nargout > 1
    Oiinvy = Yh-Pp*(Miysq);
    O = Pp*Qp';
    sigsqr = qf/Nd;
end
end

