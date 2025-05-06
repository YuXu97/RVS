function [obj, qdr, lgdt] = nglglklhd_regN(hyper, Psi, Y, kernel, M, method)
warning('off', 'MATLAB:nearlySingularMatrix');

[Nd, ~] = size(Psi);

O = CalculateOutputKernel_regN(Psi, Psi, M, kernel, hyper, 1);
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
            qdr = 1/eps;
            lgdt = 1/eps;
        else
            Linv = eye(Nd)/L;
            Oiinv =  Linv'*Linv;
            h0 = (Y'*sum(Oiinv,2))/sum(sum(Oiinv));
            Y = Y - h0;
            qdr = Nd*log(sumsqr(Linv*Y));
            lgdt = 2*sum(log(abs(diag(L))));
            obj = Nd*(1-log(Nd))+ qdr + lgdt;
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


