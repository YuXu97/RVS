function [EstInfo] = rvs_ml(data, n, M, kernel, method)

if strcmp(kernel, 'DC_mpoly')
    nstart = n;
else
    nstart = 2*n;
end

[lb,ub] = lubounds(kernel, M);
u = data(:,1);
N = length(u);
Nd = N-nstart+1;
y = data(nstart:N,2);
Psi = CalculatePsi(u, n);

hpini = ini_rvs(M, Psi, y, kernel, method);
ff = @(x)nglglklhd(x, Psi, y, kernel, M, method);
options = optimoptions('fmincon', 'Display', 'off');

if strcmp(kernel,'DC_mpoly')
    %linear constraint
    BD = zeros(M*(M-1)/2,(M+2)*(M-1)/2);
    for i = 1:M-1
        r1 = 1+i*(i-1)/2;
        r2 = i+i*(i-1)/2;
        c1 = i*(i+1)/2;
        c2 = i*(i+1)/2+i;
        BD(r1:r2,c1:c2) = toeplitz([1;zeros(r2-r1,1)],[1;-1;zeros(c2-c1-1,1)]);
    end
    
    A = [zeros(M*(M-1)/2,M+1),BD, zeros(M*(M-1)/2,M*(M+1)/2)];
    b = zeros(M*(M-1)/2,1);
    
    hp = fmincon(ff,hpini,A,b,[],[],lb,ub,[],options);
else
    hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);
end



% options = optimset('Display', 'off');
% hp = fminsearch(ff,hpini,options);


[cost, Oiinv] = nglglklhd(hp, Psi, y, kernel, M, method);

h0 = (y'*sum(Oiinv,2))/sum(Oiinv,'all');
W = Oiinv*(y-h0);
O = CalculateOutputKernel(Psi, Psi, M, kernel, hp, 1);
yhat = O*W + h0;
EstInfo.hpini = hpini;
EstInfo.hp = hp;
EstInfo.h0 = h0;
EstInfo.sigsqr = (y-h0)'*Oiinv*(y-h0)/Nd;
EstInfo.kernel = kernel;
EstInfo.method = method;
EstInfo.yhat = yhat;
EstInfo.W = W;
EstInfo.OutputKernel = O;
EstInfo.Psi = Psi;
EstInfo.cost = cost;
end

function [hpini] = ini_rvs(M, Psi, Y, kernel, method)


if strcmp(kernel, 'DC_mpoly')
    
    % hp = [c1 c2 ... cM lam1 lam21 lam22 ... lamM1 ...lamMM rho1 rho21 rho22 ... rhoM1 ...rhoMM]
    ratio = 10;
    lb = [-10*ones(M,1); 0.01*ones(M*(M+1)/2,1); 0.01*ones(M*(M+1)/2,1)];
    ub = [5*ones(M,1);  0.99*ones(M*(M+1)/2,1);    ones(M*(M+1)/2,1)];
    
    %linear constraint
    BD = zeros(M*(M-1)/2,(M+2)*(M-1)/2);
    for i = 1:M-1
        r1 = 1+i*(i-1)/2;
        r2 = i+i*(i-1)/2;
        c1 = i*(i+1)/2;
        c2 = i*(i+1)/2+i;
        BD(r1:r2,c1:c2) = toeplitz([1;zeros(r2-r1,1)],[1;-1;zeros(c2-c1-1,1)]);
    end
    
    A = [zeros(M*(M-1)/2,M+1),BD, zeros(M*(M-1)/2,M*(M+1)/2)];
    b = zeros(M*(M-1)/2,1);
    x0 = [zeros(M,1);0.05*(1:M*(M+1)/2)';0.05*(1:M*(M+1)/2)'];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',x0, 'Aineq',A,'bineq',b, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'DC_wh_bd')
    
    % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2]
    %     lb = [-10*ones(M,1);  0.01; 0.01; 0.01; 0.01];
    %     ub = [5*ones(M,1);  0.99; 0.99;  1;  1];
    %     div = [10*ones(M,1);5*ones(4,1)];
    %     fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
    %     [hpini,~,~] = GridSearch(length(lb), lb, ub, lb+(ub-lb).*0.333, div, 1e-3*ones(length(lb),1), 1e-5, 1000, fun);
    %
    
    ratio = 10;
    lb = [-10*ones(M,1);  0.01; 0.01; 0.01; 0.01];
    ub = [5*ones(M,1);  0.99; 0.99;  1;  1];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'DC_wh_opt') || strcmp(kernel, 'DC_wh_dc') ...
        || strcmp(kernel, 'DC_wh_ob') || strcmp(kernel, 'DC_wh_oba')
    
    
    % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2]
    %     lb = [-5*ones(M,1);  0.01; 0.01; 0.01; 0.01];
    %     ub = [5*ones(M,1);  0.99; 0.99;  1;  1];
    %     div = [10*ones(M,1);5*ones(4,1)];
    %     fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
    %     [hpini,~,~] = GridSearch(length(lb), lb, ub, lb+(ub-lb).*0.333, div, 1e-3*ones(length(lb),1), 1e-5, 1000, fun);
    %
    
    ratio = 10;
    lb = [-5*ones(M,1); 0.01; 0.01; 0.01; 0.01];
    ub = [5*ones(M,1);  0.99; 0.99;  1;  1];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if  strcmp(kernel, 'DC_wh_dcs')
    % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2 alpha]
    
    ratio = 10;
    lb = [-40*ones(M,1);  0.01; 0.01; 0.01; 0.01; 0.0001];
    ub = [40*ones(M,1);  0.99; 0.99;  1;  1; 0.999];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end
end

function [obj, Oiinv] = nglglklhd(hyper, Psi, Y, kernel, M, method)
warning('off', 'MATLAB:nearlySingularMatrix');

[Nd, n] = size(Psi);
if ~strcmp(kernel, 'DC_mpoly')
    Nd = Nd - n;
end

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
            h0 = (Y'*sum(Oiinv,2))/sum(Oiinv,'all');
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
            h0 = (Y'*sum(Oiinv,2))/sum(Oiinv,'all');
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
                h0 = (Y'*sum(Oiinv,2))/sum(Oiinv,'all');
                Y = Y - h0;
                obj = Nd*(1-log(Nd))+ + Nd*log(sumsqr(Linv*Y))+2*sum(log(abs(diag(L))));
                
            catch
                [U,S,~] = svd(Oi);
                Sinv = diag(1./diag(S));
                Oiinv = U*Sinv*U';
                h0 = (Y'*sum(Oiinv,2))/sum(Oiinv,'all');
                Y = Y - h0;
                obj = Nd*(1-log(Nd))+ + Nd*log(sumsqr(sqrt(Sinv)*U'*Y))+sum(log(abs(diag(S))));
                
            end
            
        end
end
end
