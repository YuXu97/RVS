function [O, V1, V2] = CalculateOutputKernel(Psi1, Psi2, M, kernel, hyper, issym)
%issym: use chol or svd to enforce the Output kernel to be symmetric (1) or
%not (0).
% in this case, Psi1 = Psi2
%O: N1 x N2 matrix
%if O is a rank-1 (or low rank) kernel matrix: O=V1*V2'
[N1, ~] = size(Psi1);
[N2, n] = size(Psi2);

O = zeros(N1,N2);
switch kernel
    case 'DC-bd-sp'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        rho = 0.99;

        ts = 1:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);

        if issym == 1
            try
                Kc = chol(K)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
        else
            Otmp = Psi1*K*Psi2';
        end

        for i = 1:M
            O = O + c(i)*Otmp.^i;
        end




    case 'DC-dc-sp'
        % hp = [c1 c2 ... cM lam rho]
        % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
        % cm = am c^m

        ts = 1:n;
        c = hyper(1:M);%.^3;
        lam = hyper(M+1);
        rho = 0.99;


        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        t = (lam*rho).^(ts);

        if issym == 1
            xi = Psi2*t';
            XI = repmat(xi,1,N2);
            try
                Kc = chol(K)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';

        else
            xi1 = Psi1*t';
            xi2 = Psi2*t';
            XI1 = repmat(xi1,1,N2);
            XI2 = repmat(xi2,1,N1);
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


        ts = 1:n;
        c = hyper(1:M);%.^3;
        lam = hyper(M+1);
        rho = 0.99;


        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);

        r = 100;
        d = 1./(((1:r)'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = (1:r)'-0.5;
        Ut = sin(v1*v2').*repmat((lam/rho).^ts',1,r);
        cf = d;

        t = 2*Ut*cf; t = t';

        if issym == 1
            xi = Psi2*t';
            XI = repmat(xi,1,N2);
            try
                Kc = chol(K)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';

        else
            xi1 = Psi1*t';
            xi2 = Psi2*t';
            XI1 = repmat(xi1,1,N2);
            XI2 = repmat(xi2,1,N1);
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