function [Ov] = CalculateOutputKernelValidation(Psiv, Psi, M, kernel, hyper)

[N1,n] = size(Psiv);
[N2,~] = size(Psi);
Ov = zeros(N1,N2);
switch kernel

    case 'TC-bd'
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        U = lam.^(1:n)'; V = ones(n,1);
        K = tril(U*V')+triu(V*U',1);
        for m = 1:M
            Ov = Ov + c(m)*(Psiv*K*Psi').^m;
        end


    case 'DC-bd'
        c = exp(hyper(1:M));
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        K = tril(U*V')+triu(V*U',1);
        for m = 1:M
            Ov = Ov + c(m)*(Psiv*K*Psi').^m;
        end

    case 'TC-tc'
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        U = lam.^(1:n)'; V = ones(n,1);
        K = tril(U*V')+triu(V*U',1);
        Otmp = Psiv*K*Psi';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        zeta = lam.^(1:n)';
        XI1 = Psiv*zeta*ones(1,N2);
        XI2 = Psi*zeta*ones(1,N1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end

    case 'DC-dc'
        c = hyper(1:M).^3;
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        K = tril(U*V')+triu(V*U',1);
        Otmp = Psiv*K*Psi';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        zeta = (lam*rho).^(1:n)';
        XI1 = Psiv*zeta*ones(1,N2);
        XI2 = Psi*zeta*ones(1,N1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end

    case 'DC-bd-sp'
        c = exp(hyper(1:M));
        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        K = tril(U*V')+triu(V*U',1);
        for m = 1:M
            Ov = Ov + c(m)*(Psiv*K*Psi').^m;
        end

    case 'DC-dc-sp'
        c = hyper(1:M);%.^3;
        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        K = tril(U*V')+triu(V*U',1);
        Otmp = Psiv*K*Psi';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        zeta = (lam*rho).^(1:n)';
        XI1 = Psiv*zeta*ones(1,N2);
        XI2 = Psi*zeta*ones(1,N1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end

    case 'DC-ob'
        c = hyper(1:M).^3;
        lam = hyper(M+1); rho = hyper(M+2);
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        K = tril(U*V')+triu(V*U',1);
        Otmp = Psiv*K*Psi';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        r = 100; ts = 1:n;
        d = 1./(((1:r)'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = (1:r)'-0.5;
        Ut = sqrt(2)*sin(v1*v2').*repmat((lam/rho).^ts',1,r);
        cf = sqrt(2)*d;

        zeta = Ut*cf;

        XI1 = Psiv*zeta*ones(1,N2);
        XI2 = Psi*zeta*ones(1,N1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end

    case 'DC-ob-sp'
        c = hyper(1:M);%.^3;
        lam = hyper(M+1); rho = 0.99;
        U = (lam*rho).^(1:n)'; V = (lam/rho).^(1:n)';
        K = tril(U*V')+triu(V*U',1);
        Otmp = Psiv*K*Psi';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        r = 100; ts = 1:n;
        d = 1./(((1:r)'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = (1:r)'-0.5;
        Ut = sqrt(2)*sin(v1*v2').*repmat((lam/rho).^ts',1,r);
        cf = sqrt(2)*d;

        zeta = Ut*cf;

        XI1 = Psiv*zeta*ones(1,N2);
        XI2 = Psi*zeta*ones(1,N1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end
end