function [obj, qdr, lgdt] = nglglklhd_sp(hyper, Pi, Rho, y, kernel, M, n, CM, Cf, indcum)
warning('off', 'MATLAB:nearlySingularMatrix');

[Nd,~] = size(Pi);

[Pp, Qp] = CalculateOutputKernelGenerators_sp(Pi, Rho, M, n, kernel, hyper, CM, Cf, indcum);

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

lgdt = logdet(M);
qdr = Nd*log(qf);

obj = Nd*(1-log(Nd)) + qdr + lgdt;

end

