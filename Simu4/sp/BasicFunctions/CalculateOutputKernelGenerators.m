function [Pp, Qp] = CalculateOutputKernelGenerators(Pi, Rho, M, n, kernel, hyper, CM, Cf, indcum)

% [N,~] = size(Pi);

switch kernel

    case 'TC-bd'
        % hp = [c1 ... CM lam]
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        U = lam.^(1:n)'; V = ones(n,1);
        K = tril(U*V')+triu(V*U',1);
        [Ut,S,~] = svd(Rho'*K*Rho);
        L = Ut*sqrt(S);
        Ub = Pi*L;
        [Pp, Qp] = PolyGenerators(M,Ub,Ub,c, CM, Cf, indcum);

    case 'TC-tc'
        % hp = [c1 ... CM lam]
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        U = lam.^(1:n)'; V = ones(n,1);
        K = tril(U*V')+triu(V*U',1);
        [Ut,S,~] = svd(Rho'*K*Rho);
        L = Ut*sqrt(S);
        Ub = Pi*L;
        zeta = lam.^(1:n)';
        psi = Pi*(Rho'*zeta);
        [Pp, Qp] = PolyFullGenerators(M,Ub,Ub,psi,c, CM, Cf, indcum);

    case 'DC-bd'
        % hp = [c1 ... CM lam rho]
        c = exp(hyper(1:M));
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        K = tril(U*V')+triu(V*U',1);
        [Ut,S,~] = svd(Rho'*K*Rho);
        L = Ut*sqrt(S);
        Ub = Pi*L;
        [Pp, Qp] = PolyGenerators(M,Ub,Ub,c, CM, Cf, indcum);

    case 'DC-dc'
        % hp = [c1 ... CM lam rho]
        c = hyper(1:M).^3;
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        K = tril(U*V')+triu(V*U',1);
        [Ut,S,~] = svd(Rho'*K*Rho);
        L = Ut*sqrt(S);
        Ub = Pi*L;
        zeta = (lam*rho).^(1:n)';
        psi = Pi*(Rho'*zeta);
        [Pp, Qp] = PolyFullGenerators(M,Ub,Ub,psi,c, CM, Cf, indcum);


    case 'DC-bd-sp'
        % hp = [c1 ... CM lam 0.99]
        c = exp(hyper(1:M));
        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        %K = tril(U*V')+triu(V*U',1);
        Krho = semi_matrix_multiplicaiton(U',V',Rho);
        [Ut,S,~] = svd(Rho'*Krho);
        %[Ut,S,~] = svd(Rho'*K*Rho);
        L = Ut*sqrt(S);
        Ub = Pi*L;
        [Pp, Qp] = PolyGenerators(M,Ub,Ub,c, CM, Cf, indcum);

    case 'DC-dc-sp'
        c = hyper(1:M);%.^3;
        
        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        %K = tril(U*V')+triu(V*U',1);
        Krho = semi_matrix_multiplicaiton(U',V',Rho);
        [Ut,S,~] = svd(Rho'*Krho);
        %[Ut,S,~] = svd(Rho'*K*Rho);
        L = Ut*sqrt(S);
        Ub = Pi*L;


        zeta = (lam*rho).^(1:n)';
        psi = Pi*(Rho'*zeta);
        [Pp, Qp] = PolyFullGenerators(M,Ub,Ub,psi,c, CM, Cf, indcum);



    case 'DC-ob'
        % hp = [c1 ... CM lam rho]
        c = hyper(1:M).^3; ts = 1:n;
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(ts)'; V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        [Ut,S,~] = svd(Rho'*K*Rho);
        L = Ut*sqrt(S);
        Ub = Pi*L;

        r = 100; rs = 1:r;
        d = 1./((rs'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = rs'-0.5;
        Ut = sin(v1*v2').*(V*ones(1,r));%repmat(V,1,r);
        cf = d;
        zeta = 2*Ut*cf;
        psi = Pi*(Rho'*zeta);
        [Pp, Qp] = PolyFullGenerators(M,Ub,Ub,psi,c, CM, Cf, indcum);

    case 'DC-ob-sp'
        % hp = [c1 ... CM lam rho]
         c = hyper(1:M);%.^3; 
         ts = 1:n;

        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(ts)'; V = (lam/rho).^(ts)';

        %K = tril(U*V')+triu(V*U',1);
        %[Ut,S,~] = svd(Rho'*K*Rho);
        Krho = semi_matrix_multiplicaiton(U',V',Rho);
        [Ut,S,~] = svd(Rho'*Krho);
        L = Ut*sqrt(S);
        Ub = Pi*L;

        r = 100; rs = 1:r;
        d = 1./((rs'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = rs'-0.5;
        Ut = sin(v1*v2').*(V*ones(1,r));%repmat(V,1,r);
        cf = d;
        zeta = 2*Ut*cf;
        
        psi = Pi*(Rho'*zeta);
        [Pp, Qp] = PolyFullGenerators(M,Ub,Ub,psi,c, CM, Cf, indcum);
end