function [Pp, Qp] = CalculateOutputKernelGenerators(Pi, Rho, M, n, kernel, hyper, CM, Cf, indcum)

[N,~] = size(Pi);

switch kernel

    case 'TC-bd'
        % hp = [c1 ... CM lam]
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        U = lam.^(1:N)'; V = ones(N,1);
        [~,p] = size(U);
        [~,r] = size(Pi);
        Gb = create_opk([U V],[Pi Rho]);
        Ub = Gb(n:N,1:p+r); Vb = Gb(n:N,p+r+1:end);
        [Pp, Qp] = PolyGenerators(M,Ub,Vb,c, CM, Cf, indcum);

    case 'TC-tc'
        % hp = [c1 ... CM lam]
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        U = lam.^(1:N)'; V = ones(N,1);
        [~,p] = size(U);
        [~,r] = size(Pi);
        Gb = create_opk([U V],[Pi Rho]);
        Ub = Gb(n:N,1:p+r); Vb = Gb(n:N,p+r+1:end);
        zeta = lam.^(1:N)';
        %psi = Pi(n:N,:)*(Rho'*zeta);
        %psi = Psi_ir*zeta;       
        psi = trilv(Pi',Rho',zeta);
        psi = psi(n:N);
        [Pp, Qp] = PolyFullGenerators(M,Ub,Vb,psi,c, CM, Cf, indcum);

    case 'DC-bd'
        % hp = [c1 ... CM lam rho]
        c = exp(hyper(1:M));
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(1:N)'; V = (lam/rho).^(1:N)';
        [~,p] = size(U);
        [~,r] = size(Pi);
        Gb = create_opk([U V],[Pi Rho]);
        Ub = Gb(n:N,1:p+r); Vb = Gb(n:N,p+r+1:end);
        [Pp, Qp] = PolyGenerators(M,Ub,Vb,c, CM, Cf, indcum);

    case 'DC-dc'
        % hp = [c1 ... CM lam rho]
        c = hyper(1:M).^3;
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(1:N)'; V = (lam/rho).^(1:N)';
        [~,p] = size(U);
        [~,r] = size(Pi);
        Gb = create_opk([U V],[Pi Rho]);
        Ub = Gb(n:N,1:p+r); Vb = Gb(n:N,p+r+1:end);
        zeta = (lam*rho).^(1:N)';
        %psi = Pi(n:N,:)*(Rho'*zeta);
        %psi = Psi_ir*zeta;       
        psi = trilv(Pi',Rho',zeta);
        psi = psi(n:N);
        [Pp, Qp] = PolyFullGenerators(M,Ub,Vb,psi,c, CM, Cf, indcum);

    case 'DC-bd-sp'
        % hp = [c1 ... CM lam 0.99]
        c = exp(hyper(1:M));
        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(1:N)'; V = (lam/rho).^(1:N)';
        [~,p] = size(U);
        [~,r] = size(Pi);
        Gb = create_opk([U V],[Pi Rho]);
        Ub = Gb(n:N,1:p+r); Vb = Gb(n:N,p+r+1:end);
        [Pp, Qp] = PolyGenerators(M,Ub,Vb,c, CM, Cf, indcum);

    case 'DC-dc-sp'
        c = hyper(1:M);%.^3;
        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(1:N)'; V = (lam/rho).^(1:N)';
        [~,p] = size(U);
        [~,r] = size(Pi);
        Gb = create_opk([U V],[Pi Rho]);
        Ub = Gb(n:N,1:p+r); Vb = Gb(n:N,p+r+1:end);

        zeta = (lam*rho).^(1:N)';

        %psi = Pi(n:N,:)*(Rho'*zeta);
        %psi = Psi_ir*zeta; 
        psi = trilv(Pi',Rho',zeta);
        psi = psi(n:N);
        %%%%
        %         nc = 200;
        %         zeta = (lam*rho).^(1:nc)';
        %         psi = Pi(n:N,:)*(Rho(1:nc,:)'*zeta);
        
        [Pp, Qp] = PolyFullGenerators(M,Ub,Vb,psi,c, CM, Cf, indcum);


    case 'DC-ob'
        % hp = [c1 ... CM lam rho]
        c = hyper(1:M).^3; ts = 1:N;
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(ts)'; V = (lam/rho).^(ts)';
        [~,p] = size(U);
        [~,r] = size(Pi);
        Gb = create_opk([U V],[Pi Rho]);
        Ub = Gb(n:N,1:p+r); Vb = Gb(n:N,p+r+1:end);


        r = 100; rs = 1:r; 
        d = 1./((rs'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = rs'-0.5;
        Ut = sin(v1*v2').*(V*ones(1,r));%repmat(V,1,r);
        cf = d;
        zeta = 2*Ut*cf;
        %psi = Pi(n:N,:)*(Rho'*zeta);
        %psi = Psi_ir*zeta;       
        psi = trilv(Pi',Rho',zeta);
        psi = psi(n:N);
        [Pp, Qp] = PolyFullGenerators(M,Ub,Vb,psi,c, CM, Cf, indcum);

    case 'DC-ob-sp'
        % hp = [c1 ... CM lam rho]
        c = hyper(1:M);%.^3;
         ts = 1:N;
        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(ts)'; V = (lam/rho).^(ts)';
        [~,p] = size(U);
        [~,r] = size(Pi);
        Gb = create_opk([U V],[Pi Rho]);
        Ub = Gb(n:N,1:p+r); Vb = Gb(n:N,p+r+1:end);

%         nc = 200;
%         tsn = ts(1:nc);
%         r = 100; rs = 1:r;
%         d = 1./((rs'-0.5).^2)/pi^2;
%         v1 = pi*rho.^(2*tsn)'; v2 = rs'-0.5;
%         Ut = sin(v1*v2').*(V(1:nc)*ones(1,r));%repmat(V,1,r);
%         cf = d;
%         zeta = 2*(Ut*cf);
%         psi = Pi(n:N,:)*(Rho(1:nc,:)'*zeta);
        %%%%%%%%%
        r = 100; rs = 1:r;
        d = 1./((rs'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = rs'-0.5;
        Ut = sin(v1*v2').*(V*ones(1,r));%repmat(V,1,r);
        cf = d;
        zeta = 2*Ut*cf;
        %psi = Pi(n:N,:)*(Rho'*zeta);
        %psi = Psi_ir*zeta;       
        psi = trilv(Pi',Rho',zeta);
        psi = psi(n:N);

        [Pp, Qp] = PolyFullGenerators(M,Ub,Vb,psi,c, CM, Cf, indcum);
end