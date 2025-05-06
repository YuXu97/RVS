function [O, V1, V2] = CalculateOutputKernel_regN(Psi1, Psi2, M, kernel, hyper, issym)
%issym: use chol or svd to enforce the Output kernel to be symmetric (1) or
%not (0).
% in this case, Psi1 = Psi2
%O: N1 x N2 matrix
%if O is a rank-1 (or low rank) kernel matrix: O=V1*V2'
[Nd1, N1] = size(Psi1);
[Nd2, N2] = size(Psi2);

O = zeros(Nd1,Nd2);
switch kernel
    case 'DC-bd-sp'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        rho = 0.99;

        if issym == 1
            ts = 1:N2;
            U = (lam*rho).^(ts)';
            V = (lam/rho).^(ts)';
            K = tril(U*V')+triu(V*U',1);
            try
                Kc = chol(K)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
        else
            K = zeros(N1,N2);
            for i = 1:N1
                for j = 1:N2
                    K(i,j) = lam^(i+j)*rho^abs(i-j);
                end
            end
            Otmp = Psi1*K*Psi2';
        end

        for i = 1:M
            O = O + c(i)*Otmp.^i;
        end




    case 'DC-dc-sp'
        % hp = [c1 c2 ... cM lam rho]
        % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
        % cm = am c^m

        
        c = hyper(1:M);%.^3;
        lam = hyper(M+1);
        rho = 0.99;

        if issym == 1
            ts = 1:N2;
            U = (lam*rho).^(ts)';
            V = (lam/rho).^(ts)';
            K = tril(U*V')+triu(V*U',1);
            t = (lam*rho).^(ts);

            xi = Psi2*t';
            XI = repmat(xi,1,Nd2);
            try
                Kc = chol(K)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
%             Otmp = Psi2*K*Psi2';

        else
            t1 = (lam*rho).^(1:N1);
            t2 = (lam*rho).^(1:N2);
            xi1 = Psi1*t1';
            xi2 = Psi2*t2';
            XI1 = repmat(xi1,1,Nd2);
            XI2 = repmat(xi2,1,Nd1);
            K = zeros(N1,N2);
            for i = 1:N1
                for j = 1:N2
                    K(i,j) = lam^(i+j)*rho^abs(i-j);
                end
            end
            Otmp = Psi1*K*Psi2';
        end

        for i = 1:M
            O = O + c(i)^2*Otmp.^i;
        end

        for i = 1:M
            for j = i+1:M
                tmp = Otmp.^i;
                if issym == 1
                    TMP = tmp.*(XI.^(j-i));
                    O = O + c(i)*c(j)*(TMP+TMP');
                else
                    O = O + c(i)*c(j)*tmp.*((XI2.^(j-i))'+XI1.^(j-i));
                end
            end
        end


    case 'DC-ob-sp'
        % hp = [c1 c2 ... cM lam rho]
        
        c = hyper(1:M);%.^3;
        lam = hyper(M+1);
        rho = 0.99;

        if issym == 1
            ts = 1:N2;
            U = (lam*rho).^(ts)';
            V = (lam/rho).^(ts)';
            K = tril(U*V')+triu(V*U',1);

            r = 100;
            d = 1./(((1:r)'-0.5).^2)/pi^2;
            v1 = pi*rho.^(2*ts)'; v2 = (1:r)'-0.5;
            Ut = sin(v1*v2').*repmat((lam/rho).^ts',1,r);
            cf = d;

            t = 2*Ut*cf; t = t';


            xi = Psi2*t';
            XI = repmat(xi,1,Nd2);
            try
                Kc = chol(K)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
%             Otmp = Psi2*K*Psi2';

        else
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

            XI1 = Psi1*zeta1*ones(1,Nd2);
            XI2 = Psi2*zeta2*ones(1,Nd1);
            K = zeros(N1,N2);
            for i = 1:N1
                for j = 1:N2
                    K(i,j) = lam^(i+j)*rho^abs(i-j);
                end
            end
            Otmp = Psi1*K*Psi2';
        end

        for i = 1:M
            O = O + c(i)^2*Otmp.^i;
        end

        for i = 1:M
            for j = i+1:M
                tmp = Otmp.^i;
                if issym == 1
                    TMP = tmp.*(XI.^(j-i));
                    O = O + c(i)*c(j)*(TMP+TMP');
                else
                    O = O + c(i)*c(j)*tmp.*((XI2.^(j-i))'+XI1.^(j-i));
                end
            end
        end

    

end




end