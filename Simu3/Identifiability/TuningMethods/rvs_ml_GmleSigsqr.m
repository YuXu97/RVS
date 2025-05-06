function [EstInfo] = rvs_ml_GmleSigsqr(data, n, M, kernel, method)

[lb,ub] = lubounds_GmleSigsqr(kernel, M);

u = data(:,1); 
N = length(u);
y = data(n:N,2);
Psi = CalculatePsi(u, n);

hpini = ini_rvs(M, Psi, y, kernel, method);
ff = @(x)nglglklhd(x, Psi, y, kernel, M, method);

%options = optimoptions('fmincon', 'Display', 'off','FunctionTolerance', 1e-7, 'StepTolerance', 1e-11);
options = optimoptions('fmincon', 'Display', 'off');
hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);

% options = optimset('Display', 'off');
% hp = fminsearch(ff,hpini,options);

[cost, Mx1, Mx2] = nglglklhd(hp, Psi, y, kernel, M, method);


% h0 = hp(end);
h0 = 0;
poly = hp(1:M); 

B = zeros(N-n+1,M);
for i = 1:M
    B(:,i) = u(n:N).^i;
end
b = B*poly;

switch method
    case 'chol'
        Oiinv =  Mx1'*Mx1;
    case 'svd'
        Oiinv =  Mx1*Mx2*Mx1';
    case 'combined'
        if Mx2(2,1) == 0
            Oiinv =  Mx1*Mx2*Mx1';
        else
            Oiinv =  Mx1'*Mx1;
        end
end




W = Oiinv*(y-h0-b);
[O, K, O1h] = CalculateOutputKernel(Psi, Psi, M, kernel, hp, 1);
yhat = O*W + h0 + b;
ghat = [1;O1h*W];

sigsqrhat = (y-h0-b)'*Oiinv*(y-h0-b)/(N-n+1);

EstInfo.b = b;
EstInfo.K = K;
EstInfo.ghat = ghat;
EstInfo.hpini = hpini;
EstInfo.hp = hp;
EstInfo.kernel = kernel;
EstInfo.method = method;
EstInfo.yhat = yhat;
%EstInfo.efit = gof(y,yhat);
EstInfo.W = W;
EstInfo.OutputKernel = O;
EstInfo.Psi = Psi;
EstInfo.cost = cost;
EstInfo.sigsqrhat = sigsqrhat;
end

function [hpini] = ini_rvs(M, Psi, y, kernel, method)

if strcmp(kernel, 'TC_bd') || strcmp(kernel, 'TC_opt+') || strcmp(kernel, 'TC_tc+')...
        || strcmp(kernel, 'TC_opt-') || strcmp(kernel, 'TC_tc-')
        
        
    % hp = [a1 a2 ... aM c lam]
    ratio = 10;
    lb = [-4*ones(M,1); -25;  0.01];
    ub = [4*ones(M,1); 10; 0.99];
    fun = @(x)nglglklhd(x, Psi, y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.4, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));

end

if strcmp(kernel, 'DC_bd') || strcmp(kernel, 'NEW_bd')...
        || strcmp(kernel, 'DC_opt+') || strcmp(kernel, 'DC_dc+') || strcmp(kernel, 'NEW_new+')...
        || strcmp(kernel, 'DC_opt-') || strcmp(kernel, 'DC_dc-') || strcmp(kernel, 'NEW_new-')...      
        || strcmp(kernel, 'DC_eig+') || strcmp(kernel, 'DC_eig-')
    % hp = [a1 a2 ... aM c lam rho]
    ratio = 10;    
    lb = [-4*ones(M,1); -25;  0.01; -1];
    ub = [4*ones(M,1); 15; 0.99; 1];
    fun = @(x)nglglklhd(x, Psi, y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));

end


end

function [obj, Mx1, Mx2] = nglglklhd(hyper, Psi, Y, kernel, M, method)
warning('off', 'MATLAB:nearlySingularMatrix');


[Nd, n] = size(Psi);


h0 = 0; %h0 = hyper(end); 
poly = hyper(1:M);

B = zeros(Nd,M);
for i = 1:M
    B(:,i) = Psi(:,1).^i;
end
b = B*poly;


Y = Y - h0 - b;
O = CalculateOutputKernel(Psi, Psi, M, kernel, hyper, 1);
Oi = O+eye(Nd);

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
            %fprintf('chol problem!\n');
            chol_flag = 0;
        end
        
        
        if  cond_flag ~= 1 || chol_flag ~= 1
            obj = 1/eps;
        else
            Linv = eye(Nd)/L;
            obj = Nd*log(sumsqr(Linv*Y))+2*sum(log(abs(diag(L)))) + Nd*(1-log(Nd));
            Mx1 = Linv; Mx2 = L;
        end
        
    case 'svd'
        if  cond_flag ~= 1
            obj = 1/eps;
        else
            [U,S,~] = svd(Oi);
            Sinv = diag(1./diag(S));
            obj = Nd*log(sumsqr(sqrt(Sinv)*U'*Y))+sum(log(abs(diag(S))))+ Nd*(1-log(Nd));
            Mx1 = U; Mx2 = Sinv;
        end
        
    case 'combined'
        
        if  cond_flag ~= 1
            obj = 1/eps;
        else
            try
                L = chol(Oi)';
                Linv = eye(Nd)/L;
                obj = Nd*log(sumsqr(Linv*Y))+2*sum(log(abs(diag(L))))+ Nd*(1-log(Nd));
                Mx1 = Linv; Mx2 = L;
            catch
                [U,S,~] = svd(Oi);
                Sinv = diag(1./diag(S));
                obj = Nd*log(sumsqr(sqrt(Sinv)*U'*Y))+sum(log(abs(diag(S))))+ Nd*(1-log(Nd));
                Mx1 = U; Mx2 = Sinv;
            end
        end
end

end

