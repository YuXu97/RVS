function [EstInfo] = rvs_ml(data, n, M, kernel, method)
%Error model: mode = {'error', wc, ph}.
%Full model: mode = {'full',[],[]}.

[lb,ub] = lubounds(kernel, M);

u = data(:,1);
N = length(u);
y = data(n:N,2);
Psi = CalculatePsi(u, n);
%
if strcmp(kernel, 'DC-bd-odd-r1')
    hpini = ini_rvs(M, Psi, y, kernel, method);
    ff = @(x)nglglklhd(x, Psi, y, kernel, M, method);
    options = optimoptions('fmincon', 'Display', 'off');
    sigsqr = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);
    
    [cost, Oiinv, U] = nglglklhd_r1(Psi, sigsqr, kernel, y, M);
    e = ones(N-n+1,1);
    Ute = U'*e;
    h0 = (sum(y)-(U'*y)'/Oiinv*Ute)/(N-n+1-Ute'/Oiinv*Ute);
    y0 = y-h0;
    W = (U'*y0-(U'*U)/(sigsqr*eye(M)+U'*U)*(U'*y0))/sigsqr;
    yhat =U*W +h0;
    EstInfo.hpini = hpini;
    EstInfo.hp = sigsqr;
    EstInfo.h0 = h0;
    EstInfo.sigsqr = sigsqr;
    EstInfo.kernel = kernel;
    EstInfo.method = method;
    EstInfo.yhat = yhat;
    EstInfo.W = W;
    EstInfo.Psi = Psi;
    EstInfo.cost = cost;
else
    
    hpini = ini_rvs(M, Psi, y, kernel, method);
    
    ff = @(x)nglglklhd(x, Psi, y, kernel, M, method);
    
    options = optimoptions('fmincon', 'Display', 'off');
    
    hp = fmincon(ff,hpini,[],[],[],[],lb,ub,[],options);
    
    
    

    
    % hp = [2.34304190766050 -2.51441753277124 1.67027195631447 ...
    %     -0.591409073750925 0.104839744000900 -0.00732245172673269 5.96471729000931 0.896892159024929 0.996202022665384];
    %
    % hp = [2.34304190766050 -2.51441753277124 1.67027195631447 -0.591409073750925...
    %     0.104839744000900 5.96471729000931 0.896892159024929 0.996202022665384];
    
    % hp = [2.34304190766050 -2.51441753277124 1.67027195631447 -0.591409073750925...
    %      0.104839744000900 3.296692456299424 0.927085753029509 0.997485311139926];
    
    %  hp = [2.34304190766050 -2.51441753277124 1.67027195631447 3.296692456299424 0.927085753029509 0.997485311139926];
    %  %
    
    % hp = [2.34304190766050 -2.51441753277124 1.67027195631447 ...
    %      -0.591409073750925 0.104839744000900 -0.00732245172673269 -0.083321191501856 0.542896916421009 0.999999985043651];
    %
%     options = optimset('Display', 'off');
%     hp = fminsearch(ff,hpini,options);
    
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
end

function [hpini] = ini_rvs(M, Psi, Y, kernel, method)

if strcmp(kernel, 'AMLS2os-bd')
    
    % hp = [c1 c2 ... cM lam rho al]
    ratio = 10;
    lb = [-5*ones(M,1);  0.01; -0.99; 0.1];
    ub = [5*ones(M,1);  0.99; 0.99; 10];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'AMLS2od-bd')
    
    % hp = [c1 c2 ... cM lam rho omg epsi]
    ratio = 10;
    lb = [-5*ones(M,1);  0.01; -0.99; 0.1; 1e-5];
    ub = [5*ones(M,1);  0.99; 0.99; 10; 1e-4];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'DC-bd')
    
    % hp = [c1 c2 ... cM lam rho].t
    ratio = 20;
    lb = [-5*ones(M,1);  0.001; 0.001];
    ub = [2*ones(M,1);  0.999; 0.999];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel', true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end




if strcmp(kernel, 'DC-bd-odd') || strcmp(kernel, 'tDC-bd-odd')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 10;
    lb = [-(M:-1:1)'; -5;  0.5; 0.05];
    ub = [(M:-1:1)'; 3; 0.9; 0.95];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
    
end

if strcmp(kernel, 'hDC-bd-odd')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 20;
    lb = [-0.5*(M:-1:1)'; -5; 0.05; 0.05; 0.05; 1e-5; 1e-5; 1e-5];
    ub = [0.5*(M:-1:1)'; 3; 0.95; 0.95; 0.95; 2; 2; 1];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
    
end

if strcmp(kernel, 'hDC-bd')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 20;
    lb = [-2*ones(M,1); -5; 0.05; 0.05; 0.05; 1e-5; 1e-5; 1e-5];
    ub = [2*ones(M,1); 3; 0.95; 0.95; 0.95; 2; 2; 1];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
    
end

if strcmp(kernel,'hDC_2dc-bd-odd-polyfix')
    ratio = 10;
    lb = [ -5; -5; 0.05; 0.05; 0.05; 0.05; 0.05; 1e-5; 1e-5; 1e-5];
    ub = [  3;  3; 0.95; 0.95; 0.95; 0.95; 0.95;    2;    2;    1];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'SI2od_dc-bd')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 25;
    lb = [-0.5*ones(M,1); -4; 0.01; 0.01; 0.01;  0.01; 0.01];
    ub = [0.5*ones(M,1); 2; 0.99; 0.99; 0.99; 0.99; 0.99];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
    
end

if strcmp(kernel, 'SI2od_dc-bd-odd')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 20;
    lb = [-0.5*(M:-1:1)'; -6; 0.05; 0.05; 0.05;  0.05; 0.05];
    ub = [0.5*(M:-1:1)'; 2; 0.95; 0.95; 0.95; 0.95; 0.95];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
    
end
if strcmp(kernel, 'SI2od_dc-bd-odd-polyfix')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 10;
    lb = [ -5;   0.001;0.001; 0.001;  0.001; 0.001];
    ub = [ 3; 0.999; 0.999; 0.999; 0.999; 0.999];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
    
end
if strcmp(kernel, 'hDC-bd-odd-polyfix')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 10;
    lb = [-5; 0.01; 0.01; 0.01; 1e-5; 1e-5; 1e-5];
    ub = [ 3; 0.99; 0.99; 0.99; 2; 2; 1];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
    
end

if strcmp(kernel, 'hOS')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 20;

    lb = [-5; 0.05; 1e-5; 1e-5; 1e-5];
    ub = [3;  0.95; 2; 2; 1];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));   
end

if strcmp(kernel, 'SI2od')
%     
    ratio = 30;
    lb = [-6; 0.5; 0.1; 0.1];
    ub = [1; 0.96; 0.95; 0.95];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));  

%         beta = 0.5:0.03:0.96;
%         betahpr = 0.1:0.05:0.95;
%         gamma = [0.1:0.1:0.9 0.95];
%         hpr = -6:1:1;
%         
%         obj = zeros(size(beta,2),size(betahpr,2),size(gamma,2),size(hpr,2));
%         for nj = 1:size(beta,2)
%             for nr = 1:size(betahpr,2)
%                 for ng = 1:size(gamma,2)
%                     for nm = 1:size(hpr,2)
%                         hpstart = [hpr(nm) beta(nj) betahpr(nr) gamma(ng)]';
%                         obj(nj,nr,ng,nm) = nglglklhd(hpstart, Psi, Y, kernel, M, method);
%                     end
%                 end
%             end
%         end
%         [~, ind_bbeta] = min(obj(:));
%         [indj, indr, indc, indd] = ind2sub([nj nr ng nm],ind_bbeta);
%         hpini = [hpr(indd) beta(indj) betahpr(indr) gamma(indc)]';
    
end

if strcmp(kernel, '2DC_SI2od') || strcmp(kernel,'2DC_SI2od-bd-odd-polyfix')
%     
    ratio = 10;
    lb = [-6; -6; 0.05; 0.05;  0.05;  0.05;  0.05; 0.05; 0.05];
    ub = [ 2;  2; 0.95; 0.95;  0.95;  0.95;  0.95; 0.95; 0.95];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));  
end

if strcmp(kernel, 'SI2od_dc') || strcmp(kernel, 'SI2od_tdc')
%     
    ratio = 20;
    lb = [-6; 0.05; 0.05; 0.05;  0.05; 0.05];
    ub = [2; 0.95; 0.95; 0.95; 0.95; 0.95];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub)); 
end

if strcmp(kernel, 'DC-bd-odd-r1')
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 5;
    lb = [1e-5];
    ub = [5];
    fun = @(x)nglglklhd_r1(Psi, x, kernel, Y, M);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
    
end

if strcmp(kernel, 'DC-opt') || strcmp(kernel, 'DC-dc') || strcmp(kernel, 'DC-dcp')
    
    
    % hp = [c1 c2 ... cM lam rho]
    ratio = 2;
    lb = [0.01; -5*ones(M-1,1);  0.01; 0.01];
    ub = [5; 5*ones(M-1,1);  0.99; 0.99];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if strcmp(kernel, 'DC-dcr')
    
    % hp = [c1 c2 ... cM lam rho alpha]
    ratio = 10;
    lb = [0.01; -25*ones(M-1,1);  0.01; -1; -0.999];
    ub = [25; 25*ones(M-1,1);  0.99; 0.99; 0.999];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end


if strcmp(kernel, 'DC-dcs')
    
    % hp = [c1 c2 ... cM lam rho alpha]
    ratio = 8;
    lb = [0.01; -25*ones(M-1,1);  0.01; -1; 0.001];
    ub = [25; 25*ones(M-1,1);  0.99; 0.99; 0.999];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',false);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if  strcmp(kernel, 'DC-ob') || strcmp(kernel, 'DC-eig')
    
    % hp = [c1 c2 ... cM lam rho alpha beta]
    ratio = 10;
    lb = [0.01; -5*ones(M-1,1);  0.01; 0.01];
    ub = [5; 5*ones(M-1,1);  0.99; 0.99];
    fun = @(x)nglglklhd(x, Psi, Y, kernel, M, method);
    ms = MultiStart('FunctionTolerance',1e-5,'UseParallel',true);
    ms.StartPointsToRun='bounds';
    ms.Display='off';
    problem = createOptimProblem('fmincon','x0',lb+(ub-lb).*0.333, 'objective',fun,'lb',lb,'ub',ub);
    hpini = run(ms,problem,ratio*length(ub));
end

if  strcmp(kernel, 'DC-ob-ex')
    
    % hp = [c1 c2 ... cM lam rho alpha beta]
    ratio = 10;
    lb = [0.01; -5*ones(M-1,1);  0.01; -1];
    ub = [5; 5*ones(M-1,1);  0.99; 0.99];
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


function [obj, Os, U] = nglglklhd_r1(Psi, sigsqr, kernel, Y, M)

[Nd,~] = size(Psi);
[~,~,U] = CalculateOutputKernel(Psi, Psi, M, kernel, [], 1);

uty = U'*Y;
qdf = (sumsqr(Y)-uty'/(sigsqr*eye(M)+U'*U)*uty)/sigsqr;
logd = log(det(eye(M)+U'*U/sigsqr))+Nd*log(sigsqr);

obj = qdf+logd;

if nargout > 1
    Os = sigsqr*eye(M)+U'*U;
end

end

