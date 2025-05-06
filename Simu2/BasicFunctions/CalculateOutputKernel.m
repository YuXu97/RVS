function [O] = CalculateOutputKernel(Psi1, Psi2, M, kernel, hyper, issym)
%issym: use chol or svd to enforce the Output kernel to be symmetric (1) or
%not (0).
% in this case, Psi1 = Psi2
%O: N1 x N2 matrix
[N1, ~] = size(Psi1);
[N2, n] = size(Psi2);
O = zeros(N1,N2);

switch kernel
    case 'DC_mpoly'
        % hp = [c1 c2 ... cM lam1 lam21 lam22 ... lamM1 ...lamMM rho1 rho21 rho22 ... rhoM1 ...rhoMM]
        % DC_mpoly: c1 (Psi Kdc1 Psi') + c2 (Psi Kdc21 Psi').*(Psi Kdc22 Psi') + ... + cM (Psi KdcM1 Psi').*(Psi KdcM2 Psi').*....*(Psi KdcMM Psi')
        for i = 1:M
            Ot = ones(N1,N2);
            for mi = 1:i
                idlam = M+i*(i-1)/2+mi;
                idrho = M+M*(M+1)/2+i*(i-1)/2+mi;
                U = (hyper(idlam)*hyper(idrho)).^(1:n)';
                V = (hyper(idlam)/hyper(idrho)).^(1:n)';
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
                
                Ot = Ot.*Otmp;
            end
            O = O + exp(hyper(i))*Ot;
        end
        
    case 'DC_wh_bd'
        % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2]
        % pci = ai^2c2c1^i, c1,c2>0
        ts = 1:n;
        pc = exp(hyper(1:M));
        lam1 = hyper(M+1); lam2 = hyper(M+2);
        rho1 = hyper(M+3); rho2 = hyper(M+4);
        
        U1 = (lam1*rho1).^(ts)';
        V1 = (lam1/rho1).^(ts)';
        K1 = tril(U1*V1')+triu(V1*U1',1);
        
        U2 = (lam2*rho2).^(ts)';
        V2 = (lam2/rho2).^(ts)';
        K2 = tril(U2*V2')+triu(V2*U2',1);
        
        if issym == 1
            try
                Kc = chol(K1)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K1);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
        else
            Otmp = Psi1*K1*Psi2';
        end
        
        for i = 1:M
            O = O + pc(i)*Otmp.^i;
        end
        
        try
            O = conv2fft(O,K2,'valid');
        catch
            O = conv2(O,K2,'valid');
        end
        O = O(2:end,2:end);
        
        
        
        
    case 'DC_wh_opt'
        % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2]
        % pci = ai sqrt(c2) c1^i, c1 in R, c2>0
        
        ts = 1:n;
        pc = hyper(1:M).^3;
        lam1 = hyper(M+1); lam2 = hyper(M+2);
        rho1 = hyper(M+3); rho2 = hyper(M+4);
        
        U1 = (lam1*rho1).^(ts)';
        V1 = (lam1/rho1).^(ts)';
        K1 = tril(U1*V1')+triu(V1*U1',1);
        
        U2 = (lam2*rho2).^(ts)';
        V2 = (lam2/rho2).^(ts)';
        K2 = tril(U2*V2')+triu(V2*U2',1);
        
        t=(lam1*rho1).^(ts);
        
        if issym == 1
            xi = Psi2*t';
            try
                Kc = chol(K1)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K1);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
            
        else
            xi1 = Psi1*t';
            xi2 = Psi2*t';
            Otmp = Psi1*K1*Psi2';
        end
        
        for i = 1:M
            O = O + pc(i)^2*Otmp.^i;
        end
        
        
        for i = 1:M
            for j = i+1:M
                if issym == 1
                    tmpp = (xi.^i)*(xi.^j)';
                    O = O + pc(i)*pc(j)*(tmpp+tmpp');
                else
                    O = O + pc(i)*pc(j)*((xi1.^i)*(xi2.^j)'+(xi1.^j)*(xi2.^i)');
                end
            end
        end
        
        try
            O = conv2fft(O,K2,'valid');
        catch
            O = conv2(O,K2,'valid');
        end
        O = O(2:end,2:end);
        
        
        
    case 'DC_wh_dc'
        % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2]
        % pci = ai sqrt(c2) c1^i, c1 in R, c2>0
        
        ts = 1:n;
        pc = hyper(1:M).^3;
        lam1 = hyper(M+1); lam2 = hyper(M+2);
        rho1 = hyper(M+3); rho2 = hyper(M+4);
        
        U1 = (lam1*rho1).^(ts)';
        V1 = (lam1/rho1).^(ts)';
        K1 = tril(U1*V1')+triu(V1*U1',1);
        
        U2 = (lam2*rho2).^(ts)';
        V2 = (lam2/rho2).^(ts)';
        K2 = tril(U2*V2')+triu(V2*U2',1);
        
        t=(lam1*rho1).^(ts);
        
        if issym == 1
            xi = Psi2*t';
            XI = repmat(xi,1,N2);
            try
                Kc = chol(K1)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K1);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
            
        else
            xi1 = Psi1*t';
            xi2 = Psi2*t';
            XI1 = repmat(xi1,1,N2);
            XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K1*Psi2';
        end
        
        for i = 1:M
            O = O + pc(i)^2*Otmp.^i;
        end
        
        for i = 1:M
            for j = i+1:M
                if issym == 1
                    TMP = (Otmp.^i).*(XI.^(j-i));
                    O = O + pc(i)*pc(j)*(TMP+TMP');
                else
                    O = O + pc(i)*pc(j)*((Otmp.^i).*(XI2.^(j-i))'+(Otmp.^i).*(XI1.^(j-i)));
                end
            end
        end
        
        
        try
            O = conv2fft(O,K2,'valid');
        catch
            O = conv2(O,K2,'valid');
        end
        O = O(2:end,2:end);
        
    case 'DC_wh_dcs'
        % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2 alpha]
        % pci = ai sqrt(c2) c1^i, c1 in R, c2>0
        
        ts = 1:n;
        pc = hyper(1:M).^3;
        lam1 = hyper(M+1); lam2 = hyper(M+2);
        rho1 = hyper(M+3); rho2 = hyper(M+4);
        alpha = hyper(M+5);
        
        U1 = (lam1*rho1).^(ts)';
        V1 = (lam1/rho1).^(ts)';
        K1 = tril(U1*V1')+triu(V1*U1',1);
        
        U2 = (lam2*rho2).^(ts)';
        V2 = (lam2/rho2).^(ts)';
        K2 = tril(U2*V2')+triu(V2*U2',1);
        
        t=(alpha*sqrt(1-rho1^2)/rho1+1-alpha)*...
            (lam1*rho1).^(ts);
        
        if issym == 1
            xi = Psi2*t';
            XI = repmat(xi,1,N2);
            try
                Kc = chol(K1)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K1);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
            
        else
            xi1 = Psi1*t';
            xi2 = Psi2*t';
            XI1 = repmat(xi1,1,N2);
            XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K1*Psi2';
        end
        
        for i = 1:M
            O = O + pc(i)^2*Otmp.^i;
        end
        
        for i = 1:M
            for j = i+1:M
                if issym == 1
                    TMP = (Otmp.^i).*(XI.^(j-i));
                    O = O + pc(i)*pc(j)*(TMP+TMP');
                else
                    O = O + pc(i)*pc(j)*((Otmp.^i).*(XI2.^(j-i))'+(Otmp.^i).*(XI1.^(j-i)));
                end
            end
        end
        
        try
            O = conv2fft(O,K2,'valid');
        catch
            O = conv2(O,K2,'valid');
        end
        O = O(2:end,2:end);
        
    case 'DC_wh_ob'
        % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2]
        % pci = ai sqrt(c2) c1^i, c1 in R, c2>0
        
        ts = 1:n;
        pc = hyper(1:M).^3;
        lam1 = hyper(M+1); lam2 = hyper(M+2);
        rho1 = hyper(M+3); rho2 = hyper(M+4);
        
        U1 = (lam1*rho1).^(ts)';
        V1 = (lam1/rho1).^(ts)';
        K1 = tril(U1*V1')+triu(V1*U1',1);
        
        U2 = (lam2*rho2).^(ts)';
        V2 = (lam2/rho2).^(ts)';
        K2 = tril(U2*V2')+triu(V2*U2',1);
        
        r = 100;
        d = 1./(((1:r)'-0.5).^2)/pi^2;
        v1 = pi*rho1.^(2*ts)'; v2 = (1:r)'-0.5;
        Ut = sqrt(2)*sin(v1*v2').*repmat((lam1/rho1).^ts',1,r);
        cf = sign(rho1)*sqrt(2)*d;
        
        t = Ut*cf; t = t';
        
        if issym == 1
            xi = Psi2*t';
            XI = repmat(xi,1,N2);
            try
                Kc = chol(K1)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K1);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
            
        else
            xi1 = Psi1*t';
            xi2 = Psi2*t';
            XI1 = repmat(xi1,1,N2);
            XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K1*Psi2';
        end
        
        for i = 1:M
            O = O + pc(i)^2*Otmp.^i;
        end
        
        for i = 1:M
            for j = i+1:M
                if issym == 1
                    TMP = (Otmp.^i).*(XI.^(j-i));
                    O = O + pc(i)*pc(j)*(TMP+TMP');
                else
                    O = O + pc(i)*pc(j)*((Otmp.^i).*(XI2.^(j-i))'+(Otmp.^i).*(XI1.^(j-i)));
                end
            end
        end
        
        
        try
            O = conv2fft(O,K2,'valid');
        catch
            O = conv2(O,K2,'valid');
        end
        O = O(2:end,2:end);
        
        
    case 'DC_wh_oba'
        % hp = [pc1 pc2 ... pcM lam1 lam2 rho1 rho2]
        % pci = ai sqrt(c2) c1^i, c1 in R, c2>0
        
        ts = 1:n;
        pc = hyper(1:M).^3;
        lam1 = hyper(M+1); lam2 = hyper(M+2);
        rho1 = hyper(M+3); rho2 = hyper(M+4);
        
        U1 = (lam1*rho1).^(ts)';
        V1 = (lam1/rho1).^(ts)';
        K1 = tril(U1*V1')+triu(V1*U1',1);
        
        U2 = (lam2*rho2).^(ts)';
        V2 = (lam2/rho2).^(ts)';
        K2 = tril(U2*V2')+triu(V2*U2',1);
        
        r = 100;
        d = 1./(((1:r)'-0.5).^2)/pi^2;
        v1 = pi*rho1.^(2*ts)'; v2 = (1:r)'-0.5;
        Ut = sqrt(2)*sin(v1*v2').*repmat((lam1/rho1).^ts',1,r);
        cf = sqrt(2)*d;
        
        t = Ut*cf; t = t';
        
        if issym == 1
            xi = Psi2*t';
            XI = repmat(xi,1,N2);
            try
                Kc = chol(K1)';
                tmp = Psi2*Kc;
            catch
                [Uk,Sk,~] = svd(K1);
                tmp = Psi2*Uk*sqrt(Sk);
            end
            Otmp = tmp*tmp';
            
        else
            xi1 = Psi1*t';
            xi2 = Psi2*t';
            XI1 = repmat(xi1,1,N2);
            XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K1*Psi2';
        end
        
        for i = 1:M
            O = O + pc(i)^2*Otmp.^i;
        end
        
        for i = 1:M
            for j = i+1:M
                if issym == 1
                    TMP = (Otmp.^i).*(XI.^(j-i));
                    O = O + pc(i)*pc(j)*(TMP+TMP');
                else
                    O = O + pc(i)*pc(j)*((Otmp.^i).*(XI2.^(j-i))'+(Otmp.^i).*(XI1.^(j-i)));
                end
            end
        end
        
        
        try
            O = conv2fft(O,K2,'valid');
        catch
            O = conv2(O,K2,'valid');
        end
        O = O(2:end,2:end);
end