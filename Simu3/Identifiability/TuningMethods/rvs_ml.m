function [EstInfo] = rvs_ml(data, n, M, kernel, method)
%Error model: mode = {'error', wc, ph}.
%Full model: mode = {'full',[],[]}.


[lb,ub] = lubounds(kernel, M);

u = data(:,1); 
N = length(u);
y = data(n:N,2);
Psi = CalculatePsi_delay(u, n);

hpini = ini_rvs(M, Psi, y, kernel, method);
ff = @(x)nglglklhd(x, Psi, y, kernel, M, method);

%options = optimoptions('fmincon', 'Display', 'off','FunctionTolerance', 1e-7, 'StepTolerance', 1e-11);
options = optimoptions('fmincon', 'Display', 'off');
hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);

% options = optimset('Display', 'off');
% hp = fminsearch(ff,hpini,options);

[cost, Mx1, Mx2] = nglglklhd(hp, Psi, y, kernel, M, method);


h0 = hp(end-1);
%h0 = 0;
poly = hp(1:M); 


sr = 1:M;
B = Psi(:,1).^(sr);
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
ghat = [0;1;O1h*W];

EstInfo.b = b;
EstInfo.h0 = h0;
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

end

function [hpini] = ini_rvs(M, Psi, Y, kernel, method)


if strcmp(kernel, 'SI2od_dc-bd')
        % hp = [c1 c2 ... cM lam rho]
    ratio = 10;
    pd = logical(mod([1:M]',2));
    lb = [-2*ones(M,1).*pd; -4; 0.01; 0.01; 0.01;  0.01; 0.01;-3;0.01];
    ub = [2*ones(M,1).*pd; 2; 0.99; 0.99; 0.99; 0.99; 0.99;3;3];

%     lb = [-2*ones(M,1); -4; 0.01; 0.01; 0.01;  0.01; 0.01;-3;0.01];
%     ub = [2*ones(M,1); 2; 0.99; 0.99; 0.99; 0.99; 0.99;3;2];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'TC_bd') || strcmp(kernel, 'TC_opt+') || strcmp(kernel, 'TC_tc+')...
        || strcmp(kernel, 'TC_opt-') || strcmp(kernel, 'TC_tc-')
        
        
    % hp = [a1 a2 ... aM c lam rho sig^2 a0]
    ratio = 10;
%     lb = [-1*ones(M,1); -15;  0.01; -4; 1];
%     ub = [1*ones(M,1); 15; 0.99; 4; 3];
    lb = [-4*ones(M,1); -25;  0.01; -6];
    ub = [4*ones(M,1); 10; 0.99; 3];
%     lb = [1.5; -3.5; -15;  0.01; -6];
%     ub = [2.5; -2.5; 15; 0.99; 3];
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
    % hp = [a1 a2 ... aM c lam rho sig^2 a0]
    ratio = 10;    
%     lb = [-1*ones(M,1); -15;  0.01; -1; -4; 1];
%     ub = [1*ones(M,1); 15; 0.99; 1; 4; 3];
    lb = [-4*ones(M,1); -25;  0.01; -1; -2];
    ub = [4*ones(M,1); 15; 0.99; 1; 8];
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


[Nd, ~] = size(Psi);

sigsqr = hyper(end); %sigsqr = 10^hyper(end-1); 
h0 = hyper(end-1); 
%h0 = 0; %
poly = hyper(1:M);

pd = logical(mod([1:M]',2));
poly = poly.*pd;

% sr = 1:2:2*M-1;
sr = 1:M;
B = Psi(:,1).^(sr);
b = B*poly;


Y = Y - h0 - b;
[O,~,~,inf_flag] = CalculateOutputKernel(Psi, Psi, M, kernel, hyper, 1);
if inf_flag
    obj = 1/eps;
    Mx1 = [];
    Mx2 = [];
    return
end
Oi = O+sigsqr*eye(Nd);

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
            obj = sumsqr(Linv*Y)+2*sum(log(abs(diag(L))));
            %obj = Nd^2*sumsqr(Linv'*Linv*Y)/sumsqr(Linv)^2;
            Mx1 = Linv; Mx2 = L;
        end
        
    case 'svd'
        if  cond_flag ~= 1
            obj = 1/eps;
        else
            [U,S,~] = svd(Oi);
            Sinv = diag(1./diag(S));
            %obj = Nd^2*sumsqr(Sinv*U'*Y)/sum(diag(Sinv))^2;
            obj = sumsqr(sqrt(Sinv)*U'*Y)+sum(log(abs(diag(S))));
            Mx1 = U; Mx2 = Sinv;
        end
        
    case 'combined'
        
        if  cond_flag ~= 1
            obj = 1/eps;
        else
            try
                L = chol(Oi)';
                Linv = eye(Nd)/L;
                %obj = Nd^2*sumsqr(Linv'*Linv*Y)/sumsqr(Linv)^2;
                obj = sumsqr(Linv*Y)+2*sum(log(abs(diag(L))));
                Mx1 = Linv; Mx2 = L;
            catch
                [U,S,~] = svd(Oi);
                Sinv = diag(1./diag(S));
                %obj = Nd^2*sumsqr(Sinv*U'*Y)/sum(diag(Sinv))^2;
                obj = sumsqr(sqrt(Sinv)*U'*Y)+sum(log(abs(diag(S))));
                Mx1 = U; Mx2 = Sinv;
            end
        end
end

end

