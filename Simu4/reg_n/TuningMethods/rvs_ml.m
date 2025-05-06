function [EstInfo] = rvs_ml(data, n, M, kernel, method)
% [lb,ub] = lubounds(kernel, M);

u = data(:,1);
N = length(u);
y = data(n:N,2);
Psi = CalculatePsi(u, n);

hpini = ini_rvs(M, Psi, y, kernel, method);

ff = @(x)nglglklhd(x, Psi, y, kernel, M, method);

% options = optimoptions('fmincon', 'Display','off');
% hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);

options = optimset('Display', 'off');%, 'MaxFunEvals', 3000, 'MaxIter', 3000, 'TolX', 1e-6,'TolFun',1e-6);
hp = fminsearch(ff,hpini,options);

[cost, Oiinv] = nglglklhd(hp, Psi, y, kernel, M, method);

h0 = (y'*sum(Oiinv,2))/sum(sum(Oiinv));
W = Oiinv*(y-h0);
O = CalculateOutputKernel(Psi, Psi, M, kernel, hp, 1);
yhat = O*W + h0;
EstInfo.hpini = hpini;
EstInfo.hp = hp;
EstInfo.h0 = h0;
EstInfo.sigsqr = (y-h0)'*Oiinv*(y-h0)/(N-n+1);
EstInfo.kernel = kernel;
EstInfo.method = method;
EstInfo.yhat = yhat;
EstInfo.W = W;
EstInfo.OutputKernel = O;
EstInfo.Psi = Psi;
EstInfo.cost = cost;
end

function [hpini] = ini_rvs(M, Psi, Y, kernel, method)
cdiv = 10; lamdiv = 10;
if strcmp(kernel, 'DC-bd-sp')
    
    % hp = [c1 c2 ... cM lam rho a0]
    
    c = linspace(-10,3,cdiv); Lc = length(c);
    lam = linspace(0.01,0.99,lamdiv); Llam = length(lam);
    obj = zeros(Lc,Lc,Lc,Llam);
    for nc1 = 1:Lc
        for nc2 = 1:Lc
            for nc3 = 1:Lc
                for nlam = 1:Llam
                    hpstart = [c(nc1);c(nc2);c(nc3);lam(nlam)];
                    obj(nc1,nc2,nc3,nlam) = nglglklhd(hpstart, Psi, Y, kernel, M, method);
                end
            end
        end
    end
    [~,indx] = min(obj(:));
    [indc1,indc2,indc3,indlam] = ind2sub([nc1 nc2 nc3 nlam],indx);
    hpini = [c(indc1);c(indc2);c(indc3);lam(indlam)];
    %%%%%%%%%%%%
%     ratio = 10;
%     lb = [-10*ones(M,1);  0.1];
%     ub = [5*ones(M,1);  0.9];
%     fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
%     ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
%     ms.StartPointsToRun='bounds';
%     ms.Display='off';
%     problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
%     hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'DC-dc-sp') || strcmp(kernel, 'DC-ob-sp')
       

    % hp = [c1 c2 ... cM lam rho]
    cp = linspace(0.001,8,cdiv); Lcp = length(cp);
    c = linspace(-8,8,cdiv); Lc = length(c);
    lam = linspace(0.01,0.99,lamdiv); Llam = length(lam);
    obj = zeros(Lcp,Lc,Lc,Llam);
    for nc1 = 1:Lcp
        for nc2 = 1:Lc
            for nc3 = 1:Lc
                for nlam = 1:Llam
                    hpstart = [cp(nc1);c(nc2);c(nc3);lam(nlam)];
                    obj(nc1,nc2,nc3,nlam) = nglglklhd(hpstart, Psi, Y, kernel, M, method);
                end
            end
        end
    end
    [~,indx] = min(obj(:));
    [indc1,indc2,indc3,indlam] = ind2sub([nc1 nc2 nc3 nlam],indx);
    hpini = [cp(indc1);c(indc2);c(indc3);lam(indlam)];
    %%%%%%%%%
%     ratio = 10;
%     lb = [0.01; -5*ones(M-1,1);  0.01; 0.01];
%     ub = [5; 5*ones(M-1,1);  0.99; 0.99];
%     fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
%     ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
%     ms.StartPointsToRun='bounds';
%     ms.Display='off';
%     problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
%     hpini = run(ms,problem,ratio*length(ub));
end
end

function [obj, Oiinv] = nglglklhd(hyper, Psi, Y, kernel, M, method)
warning('off', 'MATLAB:nearlySingularMatrix');

[Nd, ~] = size(Psi);

O = CalculateOutputKernel(Psi, Psi, M, kernel, hyper, 1);
Oi = O + eye(Nd);

cond_flag = 1; chol_flag = 1;
try
    if cond(Oi)>  1e+100
        cond_flag = 0;
    end
catch
    cond_flag = 0;
end

switch method
    case 'chol'
        try
            L = chol(Oi)';
        catch
            try
                L = chol(Oi+eps*eye(Nd))';
            catch
                try
                    L = chol(Oi+1e-4*eye(Nd))';
                    chol_flag = 0;
                catch
                    obj = 1/eps;
                    Oiinv = zeros(Nd,Nd);
                    return
                end
            end
        end


        if  cond_flag ~= 1 || chol_flag ~= 1
            obj = 1/eps;
        else
            Linv = eye(Nd)/L;
            Oiinv =  Linv'*Linv;
            h0 = (Y'*sum(Oiinv,2))/sum(sum(Oiinv));
            Y = Y - h0;
            obj = Nd*(1-log(Nd))+Nd*log(sumsqr(Linv*Y))+2*sum(log(abs(diag(L))));
        end

    case 'svd'
        if  cond_flag ~= 1
            obj = 1/eps;
        else
            [U,S,~] = svd(Oi);
            Sinv = diag(1./diag(S));
            Oiinv = U*Sinv*U';
            h0 = (Y'*sum(Oiinv,2))/sum(sum(Oiinv));
            Y = Y - h0;
            obj = Nd*(1-log(Nd))+ + Nd*log(sumsqr(sqrt(Sinv)*U'*Y))+sum(log(abs(diag(S))));
        end

    case 'combined'

        if  cond_flag ~= 1
            obj = 1/eps;
        else
            try
                L = chol(Oi)';
                Linv = eye(Nd)/L;
                Oiinv =  Linv'*Linv;
                h0 = (Y'*sum(Oiinv,2))/sum(sum(Oiinv));
                Y = Y - h0;
                obj = Nd*(1-log(Nd))+ + Nd*log(sumsqr(Linv*Y))+2*sum(log(abs(diag(L))));

            catch
                [U,S,~] = svd(Oi);
                Sinv = diag(1./diag(S));
                Oiinv = U*Sinv*U';
                h0 = (Y'*sum(Oiinv,2))/sum(sum(Oiinv));
                Y = Y - h0;
                obj = Nd*(1-log(Nd))+ + Nd*log(sumsqr(sqrt(Sinv)*U'*Y))+sum(log(abs(diag(S))));

            end

        end
end
end


