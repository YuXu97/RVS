function [O] = CalculateOutputKernel(Psi1, Psi2, M, kernel, hyper, issym)
%issym: use chol or svd to enforce the Output kernel to be symmetric (1) or
%not (0).
% in this case, Psi1 = Psi2
%O: N1 x N2 matrix
[N1, ~] = size(Psi1);
[N2, n] = size(Psi2);
N1 = N1 - n;
N2 = N2 - n;

O = zeros(N1,N2);
switch kernel
    
    case 'WH-DC'
        % hp = [c1 c2 lam1 lam2 lam3 rho1 rho2 rho3]
        ts = 1:n;
        c1 = exp(hyper(1)); c2 = exp(hyper(2));
        lam1 = hyper(3); lam2 = hyper(4); lam3 = hyper(5);
        rho1 = hyper(6); rho2 = hyper(7); rho3 = hyper(8);
        
        U1 = (lam1*rho1).^(ts)';
        V1 = (lam1/rho1).^(ts)';
        K1 = tril(U1*V1')+triu(V1*U1',1);
        
        U2 = (lam2*rho2).^(ts)';
        V2 = (lam2/rho2).^(ts)';
        K2 = tril(U2*V2')+triu(V2*U2',1);
        
        U3 = (lam3*rho3).^(ts)';
        V3 = (lam3/rho3).^(ts)';
        K3 = tril(U3*V3')+triu(V3*U3',1);
        
        %         if issym == 1
        %             try
        %                 Kc1= chol(K1)';
        %                 tmp1 = Psi2*Kc1;
        %                 Kc2= chol(K2)';
        %                 tmp2 = Psi2*Kc2;
        %             catch
        %                 [Uk1,Sk1,~] = svd(K1);
        %                 tmp1 = Psi2*Uk1*sqrt(Sk1);
        %                 [Uk2,Sk2,~] = svd(K2);
        %                 tmp2 = Psi2*Uk2*sqrt(Sk2);
        %             end
        %             O1 = tmp1*tmp1';
        %             Otmp = tmp2*tmp2';
        %         else
        %             O1 = Psi1*K1*Psi2';
        %             Otmp = Psi1*K2*Psi2';
        %         end
        
        O1 = Psi1(n+1:end,:)*K1*Psi2(n+1:end,:)';
        Otmp = Psi1*K2*Psi2';
        
        try
            O2 = conv2fft((Otmp).^2, K3,'valid');
        catch
            O2 = conv2((Otmp).^2, K3,'valid');
        end
        O2 = O2(2:end,2:end);
        
        
        O = c1*O1 + c2*O2;
        
end