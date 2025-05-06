function [O] = CalculateOutputKernel_NoHadamard(Psi1, Psi2, Phi2f, Phi2b, indx, M, kernel, hyper, issym)
%issym: use chol or svd to enforce the Output kernel to be symmetric (1) or
%not (0).
% in this case, Psi1 = Psi2
%O: N1 x N2 matrix
[N1, ~] = size(Psi1);
[N2, n] = size(Psi2);

O = zeros(N1,N2);
switch kernel
        
    case 'DC2-DC'
        % hp = [c1 c2 lam1 lam2 lam3 rho1 rho2 rho3]
        ts = 1:n;
        c1 = exp(hyper(1)); c2 = exp(hyper(2));
        lam1 = hyper(3); lam2 = hyper(4); lam3 = hyper(5);
        rho1 = hyper(6); rho2 = hyper(7); rho3 = hyper(8);
        
        U1 = (lam1*rho1).^(ts)';
        V1 = (lam1/rho1).^(ts)';
        K1 = tril(U1*V1')+triu(V1*U1',1);     
        O1 = Psi1*K1*Psi2';
        
        ut2 = (lam3*rho3).^(indx(:,2));
        vt2 = (lam3/rho3).^(indx(:,2));
        Pk2 = tril(ut2*vt2')+triu(vt2*ut2',1);
        
        v = lam2.^(indx(:,1));
        Pk11 = v*v';
        
        Pk12 = zeros(n*(n+1)/2,n*(n+1)/2);
        for i = 1:n*(n+1)/2
            for j = i:n*(n+1)/2
                Pk12(i,j) = rho2^abs(indx(i,1)-indx(j,1));
                Pk12(j,i) = Pk12(i,j);
            end
        end
        P = Pk11.*Pk12.*Pk2;
        
%         P = zeros(n*(n+1)/2,n*(n+1)/2);
%         for i = 1:n*(n+1)/2
%             for j = i:n*(n+1)/2
%                 t = indx(i,:);
%                 s = indx(j,:);
%                 p1 = lam2^(abs(t(1))+abs(s(1)))*rho2^abs(abs(t(1))-abs(s(1)));
%                 p2 = lam3^(abs(t(2))+abs(s(2)))*rho3^abs(abs(t(2))-abs(s(2)));
%                 P(i,j) = p1*p2;
%                 P(j,i) = P(i,j);
%             end
%         end
        if issym == 0
            O2 = Phi2f*P*Phi2b';
        else
            O2 = Phi2b*P*Phi2b';
        end
        
        O = c1*O1 + c2*O2;
end