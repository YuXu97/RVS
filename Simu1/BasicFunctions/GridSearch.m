function [XBest,BestF,Iters]=GridSearch(N, XLo, XHi, Xini, NumDiv, MinDeltaX, Eps_Fx, MaxIter, myFx)
% Function performs multivariate optimization using the
% grid search.
%
% Input
%
% N - number of variables
% XLo - array of lower values
% XHi - array of higher values
% NumDiv - array of number of divisions along each dimension
% MinDeltaX - array of minimum search values for each variable
% Eps_Fx - tolerance for difference in successive function values
% MaxIter - maximum number of iterations
% myFx - name of the optimized function
%
% Output
%
% XBest - array of optimized variables
% BestF - function value at optimum
% Iters - number of iterations
%
% EXAMPLE:
% lb = [-20*ones(M,1);  0.1];
% ub = [15*ones(M,1);  0.9];
% div = [10*ones(M,1);5];
% fun = @(x)nglglklhd(x, Pi, Rho, Y, kernel, M, n, CM, Cf, indcum);
% [hpini,~,~] = GridSearch(length(lb), lb, ub, lb+(ub-lb).*0.333, div, 1e-3*ones(length(lb),1), 1e-5, 1000, fun);
% 

Xcenter = (XHi + XLo) / 2;
XBest = Xcenter;
DeltaX = (XHi - XLo) ./ NumDiv;
% BestF = feval(myFx, XBest, N);
BestF = myFx(XBest);

if BestF >= 0
    LastBestF = BestF + 100;
else
    LastBestF = 100 - BestF;
end
X = Xini; % initial search value

Iters = 0;
bGoOn = 1;

while (bGoOn > 0) && (abs(BestF - LastBestF) > Eps_Fx) && (Iters <= MaxIter)

    bGoOn2 = 1;

    while bGoOn2 > 0

        Iters = Iters + 1;

        %     F = feval(myFx, X, N);
        F = myFx(X);
        if F < BestF
            LastBestF = BestF;
            BestF = F;
            XBest = X;
        end

        %*****************************************************
        % The next For loop implements a programming tricks
        % that simulated nested loops using just one For loop
        %*****************************************************
        % search next grid node
        for i = 1:N
            if X(i) >= XHi(i)
                if i < N
                    X(i) = XLo(i);
                else
                    bGoOn2 = 0;
                    break
                end
            else
                X(i) = X(i) + DeltaX(i);
                break
            end
        end

    end % while bGoOn2 > 0


    XCenter = XBest;
    DeltaX = DeltaX ./ NumDiv;
    XLo = XCenter - DeltaX .* NumDiv / 2;
    XHi = XCenter + DeltaX .* NumDiv / 2;
    X = XLo; % set initial X

    bGoOn = 0;
    for i=1:N
        if DeltaX(i) > MinDeltaX(i)
            bGoOn = 1;
        end
    end

end % while bGoOn > 0 && () && ()


