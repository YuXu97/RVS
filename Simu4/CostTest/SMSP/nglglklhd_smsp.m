function [obj, qdr, lgdt] = nglglklhd_smsp(hyper, Pi, Rho, Psi_ir, y, kernel, M, n, CM, Cf, indcum)
warning('off', 'MATLAB:nearlySingularMatrix');

[N,~] = size(Pi); Nd = N-n+1;

[Pp, Qp] = CalculateOutputKernelGenerators_smsp(Pi, Rho, Psi_ir, M,  n, kernel, hyper, CM, Cf, indcum);

% bUs = Pp;
% bVs = Qp;
[Q,R] = qr(Pp,0);
bUs = Q;
bVs = Qp*R';
s = sqrt(norm(bUs,'fro')/norm(bVs,'fro'));
if isinf(s)
    obj = 1/eps;
    qdr = 1/eps;
    lgdt = 1/eps;
    return
else
    bUs = bUs/s;
    bVs = bVs*s;
end

try
    [Wt,c] = egrss_potrf(bUs',bVs',1);
catch
    obj = 1/eps;
    qdr = 1/eps;
    lgdt = 1/eps;
    return
end
lgdt = 2*sum(reallog(c));

l1 = egrss_trsv(bUs',Wt,c,y);
le = egrss_trsv(bUs',Wt,c,ones(Nd,1));
h0 = l1'*le/sumsqr(le);

y = y - h0;

Liy = egrss_trsv(bUs',Wt,c,y);

qdr = Nd*log(sumsqr(Liy));

obj = Nd - Nd*log(Nd) + qdr + lgdt;

end


function [c,ceq] = nonlcon_dc(x)
c = x(end-1) - abs(x(end));
ceq = [];
end