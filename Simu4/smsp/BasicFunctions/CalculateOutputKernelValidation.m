function [Ov] = CalculateOutputKernelValidation(Psiv_ir, Psi_ir, M, kernel, hyper)

[Nd1,N1] = size(Psiv_ir);
[Nd2,N2] = size(Psi_ir);
Ov = zeros(Nd1,Nd2);
switch kernel

    case 'TC-bd'
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        K =zeros(N1,N2);
        for i = 1:N1
            for j = 1:N2
                K(i,j) = lam^max(i,j);
            end
        end
        for m = 1:M
            Ov = Ov + c(m)*(Psiv_ir*K*Psi_ir').^m;
        end


    case 'DC-bd'
        c = exp(hyper(1:M));
        lam = hyper(M+1); rho = hyper(M+2);
        K =zeros(N1,N2);
        for i = 1:N1
            for j = 1:N2
                K(i,j) = lam^(i+j)*rho^abs(i-j);
            end
        end
        for m = 1:M
            Ov = Ov + c(m)*(Psiv_ir*K*Psi_ir').^m;
        end

    case 'TC-tc'
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        K =zeros(N1,N2);
        for i = 1:N1
            for j = 1:N2
                K(i,j) = lam^max(i,j);
            end
        end
        Otmp = Psiv_ir*K*Psi_ir';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        zeta1 = lam.^(1:N1)';
        zeta2 = lam.^(1:N2)';
        XI1 = Psiv_ir*zeta1*ones(1,Nd2);
        XI2 = Psi_ir*zeta2*ones(1,Nd1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end

    case 'DC-dc'
        c = hyper(1:M).^3;
        lam = hyper(M+1); rho = hyper(M+2);
        K =zeros(N1,N2);
        for i = 1:N1
            for j = 1:N2
                K(i,j) = lam^(i+j)*rho^abs(i-j);
            end
        end
        Otmp = Psiv_ir*K*Psi_ir';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        zeta1 = (lam*rho).^(1:N1)';
        zeta2 = (lam*rho).^(1:N2)';
        XI1 = Psiv_ir*zeta1*ones(1,Nd2);
        XI2 = Psi_ir*zeta2*ones(1,Nd1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end

    case 'DC-bd-sp'
        c = exp(hyper(1:M));
        lam = hyper(M+1); rho = 0.99;
        K =zeros(N1,N2);
        for i = 1:N1
            for j = 1:N2
                K(i,j) = lam^(i+j)*rho^abs(i-j);
            end
        end
        for m = 1:M
            Ov = Ov + c(m)*(Psiv_ir*K*Psi_ir').^m;
        end

    case 'DC-dc-sp'
        c = hyper(1:M);%.^3;
        lam = hyper(M+1); rho = 0.99;
        K =zeros(N1,N2);
        for i = 1:N1
            for j = 1:N2
                K(i,j) = lam^(i+j)*rho^abs(i-j);
            end
        end
        Otmp = Psiv_ir*K*Psi_ir';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        zeta1 = (lam*rho).^(1:N1)';
        zeta2 = (lam*rho).^(1:N2)';
        XI1 = Psiv_ir*zeta1*ones(1,Nd2);
        XI2 = Psi_ir*zeta2*ones(1,Nd1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end

    case 'DC-ob'
        c = hyper(1:M).^3;
        lam = hyper(M+1); rho = hyper(M+2);
        K =zeros(N1,N2);
        for i = 1:N1
            for j = 1:N2
                K(i,j) = lam^(i+j)*rho^abs(i-j);
            end
        end
        Otmp = Psiv_ir*K*Psi_ir';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        r = 100;
        ts1 = 1:N1; ts2 = 1:N2;
        d = 1./(((1:r)'-0.5).^2)/pi^2;v2 = (1:r)'-0.5;
        v11 = pi*rho.^(2*ts1)';
        v12 = pi*rho.^(2*ts2)';

        Ut1 = sqrt(2)*sin(v11*v2').*repmat((lam/rho).^ts1',1,r);
        Ut2 = sqrt(2)*sin(v12*v2').*repmat((lam/rho).^ts2',1,r);
        cf = sqrt(2)*d;
        zeta1 = Ut1*cf;
        zeta2 = Ut2*cf;

        XI1 = Psiv_ir*zeta1*ones(1,Nd2);
        XI2 = Psi_ir*zeta2*ones(1,Nd1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end
        
    case 'DC-ob-sp'
        c = hyper(1:M);%.^3;
        lam = hyper(M+1); rho = 0.99;
        K =zeros(N1,N2);
        for i = 1:N1
            for j = 1:N2
                K(i,j) = lam^(i+j)*rho^abs(i-j);
            end
        end
        Otmp = Psiv_ir*K*Psi_ir';
        for m = 1:M
            Ov = Ov + c(m)^2*Otmp.^m;
        end

        r = 100;
        ts1 = 1:N1; ts2 = 1:N2;
        d = 1./(((1:r)'-0.5).^2)/pi^2;v2 = (1:r)'-0.5;
        v11 = pi*rho.^(2*ts1)';
        v12 = pi*rho.^(2*ts2)';

        Ut1 = sqrt(2)*sin(v11*v2').*repmat((lam/rho).^ts1',1,r);
        Ut2 = sqrt(2)*sin(v12*v2').*repmat((lam/rho).^ts2',1,r);
        cf = sqrt(2)*d;
        zeta1 = Ut1*cf;
        zeta2 = Ut2*cf;

        XI1 = Psiv_ir*zeta1*ones(1,Nd2);
        XI2 = Psi_ir*zeta2*ones(1,Nd1);
        for i = 1:M
            for j = i+1:M
                Ov = Ov + c(i)*c(j)*(Otmp.^i).*((XI2.^(j-i))'+(XI1.^(j-i)));
            end
        end
end