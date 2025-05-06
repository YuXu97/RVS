function [EstInfo] = rvs_ml_NoHadamard(data, n, M, kernel, method)

[lb,ub] = lubounds(kernel, M);

u = data(:,1);
N = length(u);
y = data(n:N,2);
Psi = CalculatePsi(u, n);

Phi21 = []; Phi22 = [];
for i = 0:n-1
    Phi21 = [Phi21 Psi(:,1:n-i)];
    if i ~= 0
        Phi22 = [Phi22 2*Psi(:,i+1:n)];
    else
        Phi22 = [Phi22 Psi(:,i+1:n)];
    end
end
Phi2 = Phi21.*Phi22;


[TI,TJ] = meshgrid(1:n, 1:n);
indxi = []; indxj = [];
for i = 0:n-1
    indxi = [indxi; diag(TI,i)];
end
for i = 0:n-1
    indxj = [indxj; diag(TJ,i)];
end
indx = [indxj,indxi];
indx = indx*[cos(pi/4), -sin(pi/4);sin(pi/4), cos(pi/4)];

hpini = ini_rvs(M, Psi, Phi2, indx, y, kernel, method);
ff = @(x)nglglklhd(x, Psi, Phi2, indx, y, kernel, M, method);
options = optimoptions('fmincon', 'Display', 'off');

hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);

% options = optimset('Display', 'off');
% hp = fminsearch(ff,hpini,options);

[cost, Oiinv] = nglglklhd(hp, Psi, Phi2, indx, y, kernel, M, method);

h0 = (y'*sum(Oiinv,2))/sum(sum(Oiinv));
W = Oiinv*(y-h0);
O = CalculateOutputKernel_NoHadamard(Psi, Psi, Phi2, Phi2, indx, M, kernel, hp, 1);
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

function [hpini] = ini_rvs(M, Psi, Phi2, indx, Y, kernel, method)

if strcmp(kernel, 'DC2-DC')
       
    % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2]
    
    ratio = 10;
    lb = [-20*ones(2,1);  0.01; 0.01; 0.01; 0.01; 0.01; 0.01];
    ub = [15*ones(2,1);  0.99; 0.99; 0.99;  1;  1; 1];
    fun = @(x)nglglklhd(x, Psi, Phi2, indx, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

end

function [obj, Oiinv] = nglglklhd(hyper, Psi, Phi2, indx, Y, kernel, M, method)
warning('off', 'MATLAB:nearlySingularMatrix');

[Nd, ~] = size(Psi);

O = CalculateOutputKernel_NoHadamard(Psi, Psi, Phi2, Phi2, indx, M, kernel, hyper, 1);
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
