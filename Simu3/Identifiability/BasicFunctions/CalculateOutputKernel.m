function [O, K, O1h, inf_flag] = CalculateOutputKernel(Psi1, Psi2, M, kernel, hyper, issym)
%issym: use chol or svd to enforce the Output kernel to be symmetric (1) or
%not (0).
% in this case, Psi1 = Psi2
%O: N1 x N2 matrix
[N1, ~] = size(Psi1);
[N2, n] = size(Psi2);
b1 = Psi1(:,1);
b2 = Psi2(:,1);

B1 = b1.^(0:M-1);
B2 = b2.^(0:M-1);







Psi1 = Psi1(:,2:n);
Psi2 = Psi2(:,2:n);

O = zeros(N1,N2);

inf_flag = 0;
switch kernel
    case 'SI2od_dc-bd'
        ts = (2:n)';
        pl = hyper(1:M);
        c = exp(hyper(M+1));

        lam = hyper(M+2);
        rho = hyper(M+3);
        
        lama = hyper(M+4);
        xi = hyper(M+5);
        lamg = hyper(M+6);
        
        al = -log(lama);
        bt = al/xi*sqrt(1-xi^2);
        gm = -log(lamg);
        amg = al-gm;
        
        if amg == 0
            inf_flag = 1;
            O = [];
            K = [];
            O1h = [];
            return
        end
        
        phi = acos(amg/sqrt(bt^2+amg^2));
        
        rl = lama.^ts;
        T = (rl*rl');
        
        rl1 = [cos(bt*ts) sqrt((al^2+1)/bt^2)*sin(bt*ts)];
        
        T1 = rl1*rl1';
        
        T2 = al*sin(bt*(ts+ts'))/bt;
        
        TE = exp(2*amg*min(ts,ts'));
        T3 = cos(bt*(ts-ts')).*(TE-1)/(4*bt^2*amg);
        
        TP = phi+bt*(ts+ts');
        T4 = cos(TP)-TE.*cos(2*bt*min(ts,ts')-TP);
        T4 = T4/(4*bt^2)/sqrt(bt^2+amg^2);
        
        K = T.*(T1+T2+T3+T4);
        
        rld = lam.^ts;
        Kdc = (rld*rld').*(rho.^abs(ts'-ts));
           
        Otmp = Psi1*(c*Kdc+K)*Psi2';
        
%         for i = 1:M
%             O = O + pl(i)^2*Otmp.^i;
%         end
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*pl(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(Otmp.^m);
        end
        
        if nargout > 2
            K = c*Kdc+K;
            O1h = K*(Psi2.*repmat(B2*((1:M)'.*pl),1,n-1))';
        end
    
    case 'TC_bd'
        % hp = [a1 a2 ... aM c lam sig^2 a0]
        % DC_bd: c1 (Psi Ktc Psi') + c2 (Psi Ktc Psi').^2 + ... + cM (Psi Ktc Psi').^M
        
        c = exp(hyper(M+1));
        
        lam = hyper(M+2);
        poly = hyper(1:M);
        
        ts = 2:n;
        t=lam.^ts;
        triuP = triu(repmat(t,n-1,1),1);
        K = triuP + triuP' + diag(t);
        
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
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        if nargout > 2
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
        end
        
        
    case 'TC_opt+'
        %hp = [a1 a2 ... aM c lam sig^2 a0]
        % TC_full: c1 (Psi Ktc Psi') + c2 (Psi Ktc Psi').^2 + ... +
        % cM (Psi Ktc Psi').^M + off-diagonal blocks
        
        c = exp(hyper(M+1)); lam = hyper(M+2);
        poly = hyper(1:M);
        
        ts = 2:n;
        t=lam.^ts;
        triuP = triu(repmat(t,n-1,1),1);
        K = triuP + triuP' + diag(t);
        
        
        
        xi1 = Psi1*t';
        xi2 = Psi2*t';
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
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        for i = 1:M
            for j = i+1:M
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    tmp1i = (xi1.^i).*(bt1*cc1);
                    tmp2j = (xi2.^j).*(bt2*cc2);
                    TMP = tmp1i*tmp2j';
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    tmp1i = (xi1.^i).*(bt1*cc1);
                    tmp2j = (xi2.^j).*(bt2*cc2);
                    tmp1j = (xi1.^j).*(bt1j*cc2);
                    tmp2i = (xi2.^i).*(bt2i*cc1);
                    O = O + c^(i+j)*(tmp1i*tmp2j'+tmp1j*tmp2i');
                end
                
                
            end
        end
        
        if nargout > 2
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*t'*((bt2*cc2).*xi2.^j)';
            end
            
        end
        
    case 'TC_opt-'
        %hp = [a1 a2 ... aM c lam sig^2 a0]
        % TC_full: c1 (Psi Ktc Psi') + c2 (Psi Ktc Psi').^2 + ... +
        % cM (Psi Ktc Psi').^M + off-diagonal blocks
        
        c = -exp(hyper(M+1)); lam = hyper(M+2);
        poly = hyper(1:M);
        
        ts = 2:n;
        t=lam.^ts;
        triuP = triu(repmat(t,n-1,1),1);
        K = triuP + triuP' + diag(t);
        
        
        
        xi1 = Psi1*t';
        xi2 = Psi2*t';
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
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        for i = 1:M
            for j = i+1:M
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    tmp1i = (xi1.^i).*(bt1*cc1);
                    tmp2j = (xi2.^j).*(bt2*cc2);
                    TMP = tmp1i*tmp2j';
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    tmp1i = (xi1.^i).*(bt1*cc1);
                    tmp2j = (xi2.^j).*(bt2*cc2);
                    tmp1j = (xi1.^j).*(bt1j*cc2);
                    tmp2i = (xi2.^i).*(bt2i*cc1);
                    O = O + c^(i+j)*(tmp1i*tmp2j'+tmp1j*tmp2i');
                end
                
                
            end
        end
        
        if nargout > 2
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*t'*((bt2*cc2).*xi2.^j)';
            end
            
        end
        
        
    case 'TC_tc+'
        % hp = [a1 a2 ... aM c lam sig^2 a0]
        % TC_full: c1 (Psi Ktc Psi') + c2 (Psi Ktc Psi').^2 + ... +
        % cM (Psi Ktc Psi').^M + off-diagonal blocks
        
        c = exp(hyper(M+1)); lam = hyper(M+2);
        poly = hyper(1:M);
        
        ts = 2:n;
        t=lam.^ts;
        triuP = triu(repmat(t,n-1,1),1);
        K = triuP + triuP' + diag(t);
        
        
        if issym == 1
            xi = Psi2*t';
            %             XI = repmat(xi,1,N2);
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
            %             XI1 = repmat(xi1,1,N2);
            %             XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K*Psi2';
        end
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        
        
        for i = 1:M
            for j = i+1:M
                
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    TMP = ((bt1*cc1)*((bt2*cc2).*xi.^(j-i))').*(Otmp.^i);
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    O = O + c^(i+j)*(((bt1*cc1)*((bt2*cc2).*xi2.^(j-i))').*(Otmp.^i)+...
                        (((bt1j*cc2).*xi1.^(j-i))*(bt2i*cc1)').*(Otmp.^i));
                end
            end
        end
        
        if nargout > 2
            
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*K*(Psi2.*repmat((bt2*cc2).*xi.^(j-1),1,n-1))';
            end
            
        end
        
    case 'TC_tc-'
        % hp = [a1 a2 ... aM c lam sig^2 a0]
        % TC_full: c1 (Psi Ktc Psi') + c2 (Psi Ktc Psi').^2 + ... +
        % cM (Psi Ktc Psi').^M + off-diagonal blocks
        
        c = -exp(hyper(M+1)); lam = hyper(M+2);
        poly = hyper(1:M);
        
        ts = 2:n;
        t=lam.^ts;
        triuP = triu(repmat(t,n-1,1),1);
        K = triuP + triuP' + diag(t);
        
        
        if issym == 1
            xi = Psi2*t';
            %             XI = repmat(xi,1,N2);
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
            %             XI1 = repmat(xi1,1,N2);
            %             XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K*Psi2';
        end
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        
        
        for i = 1:M
            for j = i+1:M
                
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    TMP = ((bt1*cc1)*((bt2*cc2).*xi.^(j-i))').*(Otmp.^i);
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    O = O + c^(i+j)*(((bt1*cc1)*((bt2*cc2).*xi2.^(j-i))').*(Otmp.^i)+...
                        (((bt1j*cc2).*xi1.^(j-i))*(bt2i*cc1)').*(Otmp.^i));
                end
            end
        end
        
        if nargout > 2
            
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*K*(Psi2.*repmat((bt2*cc2).*xi.^(j-1),1,n-1))';
            end
            
        end
        
        
        
    case 'DC_bd'
        % hp = [a1 a2 ... aM c lam rho sig^2 a0]
        % DC_bd: c1 (Psi Kdc Psi') + c2 (Psi Kdc Psi').^2 + ... + cM (Psi Kdc Psi').^M
        c = exp(hyper(M+1)); lam = hyper(M+2); rho = hyper(M+3);
        poly = hyper(1:M);
        
        ts = 2:n;
        
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
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        if nargout > 2
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
        end
        
        
    case 'DC_opt+'
        %hp = [a1 a2 ... aM c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi Kdc Psi') + c2 (Psi Kdc Psi').^2 + ... +
        % cM (Psi Kdc Psi').^M + off-diagonal blocks
        
        c = exp(hyper(M+1)); lam = hyper(M+2); rho = hyper(M+3);
        poly = hyper(1:M);
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        
        t=(lam*rho).^(ts);
        
        xi1 = Psi1*t';
        xi2 = Psi2*t';
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
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        for i = 1:M
            for j = i+1:M
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    tmp1i = (xi1.^i).*(bt1*cc1);
                    tmp2j = (xi2.^j).*(bt2*cc2);
                    TMP = tmp1i*tmp2j';
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    tmp1i = (xi1.^i).*(bt1*cc1);
                    tmp2j = (xi2.^j).*(bt2*cc2);
                    tmp1j = (xi1.^j).*(bt1j*cc2);
                    tmp2i = (xi2.^i).*(bt2i*cc1);
                    O = O + c^(i+j)*(tmp1i*tmp2j'+tmp1j*tmp2i');
                end
                
                
            end
        end
        
        if nargout > 2
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*t'*((bt2*cc2).*xi2.^j)';
            end
            
        end
        
    case 'DC_opt-'
        %hp = [a1 a2 ... aM c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi Kdc Psi') + c2 (Psi Kdc Psi').^2 + ... +
        % cM (Psi Kdc Psi').^M + off-diagonal blocks
        
        c = -exp(hyper(M+1)); lam = hyper(M+2); rho = hyper(M+3);
        poly = hyper(1:M);
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        
        t=(lam*rho).^(ts);
        
        xi1 = Psi1*t';
        xi2 = Psi2*t';
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
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        for i = 1:M
            for j = i+1:M
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    tmp1i = (xi1.^i).*(bt1*cc1);
                    tmp2j = (xi2.^j).*(bt2*cc2);
                    TMP = tmp1i*tmp2j';
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    tmp1i = (xi1.^i).*(bt1*cc1);
                    tmp2j = (xi2.^j).*(bt2*cc2);
                    tmp1j = (xi1.^j).*(bt1j*cc2);
                    tmp2i = (xi2.^i).*(bt2i*cc1);
                    O = O + c^(i+j)*(tmp1i*tmp2j'+tmp1j*tmp2i');
                end
                
                
            end
        end
        
        if nargout > 2
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*t'*((bt2*cc2).*xi2.^j)';
            end
            
        end
        
        
        
    case 'DC_dc+'
        % hp = [a1 a2 ... aM c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi Kdc Psi') + c2 (Psi Kdc Psi').^2 + ... +
        % cM (Psi Kdc Psi').^M + off-diagonal blocks
        
        c = exp(hyper(M+1)); lam = hyper(M+2); rho = hyper(M+3);
        poly = hyper(1:M);
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        
        t=(lam*rho).^(ts);
        
        if issym == 1
            xi = Psi2*t';
            %             XI = repmat(xi,1,N2);
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
            %             XI1 = repmat(xi1,1,N2);
            %             XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K*Psi2';
        end
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        
        
        for i = 1:M
            for j = i+1:M
                
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    TMP = ((bt1*cc1)*((bt2*cc2).*xi.^(j-i))').*(Otmp.^i);
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    O = O + c^(i+j)*(((bt1*cc1)*((bt2*cc2).*xi2.^(j-i))').*(Otmp.^i)+...
                        (((bt1j*cc2).*xi1.^(j-i))*(bt2i*cc1)').*(Otmp.^i));
                end
            end
        end
        
        if nargout > 2
            
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*K*(Psi2.*repmat((bt2*cc2).*xi.^(j-1),1,n-1))';
            end
            
        end
        
    case 'DC_dc-'
        % hp = [a1 a2 ... aM c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi Kdc Psi') + c2 (Psi Kdc Psi').^2 + ... +
        % cM (Psi Kdc Psi').^M + off-diagonal blocks
        
        c = -exp(hyper(M+1)); lam = hyper(M+2); rho = hyper(M+3);
        poly = hyper(1:M);
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        
        t=(lam*rho).^(ts);
        
        if issym == 1
            xi = Psi2*t';
            %             XI = repmat(xi,1,N2);
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
            %             XI1 = repmat(xi1,1,N2);
            %             XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K*Psi2';
        end
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        
        
        for i = 1:M
            for j = i+1:M
                
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    TMP = ((bt1*cc1)*((bt2*cc2).*xi.^(j-i))').*(Otmp.^i);
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    O = O + c^(i+j)*(((bt1*cc1)*((bt2*cc2).*xi2.^(j-i))').*(Otmp.^i)+...
                        (((bt1j*cc2).*xi1.^(j-i))*(bt2i*cc1)').*(Otmp.^i));
                end
            end
        end
        
        if nargout > 2
            
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*K*(Psi2.*repmat((bt2*cc2).*xi.^(j-1),1,n-1))';
            end
            
        end
        
    case 'DC_eig+'
        % hp = [a1 a2 ... aM c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi Kdc Psi') + c2 (Psi Kdc Psi').^2 + ... +
        % cM (Psi Kdc Psi').^M + off-diagonal blocks
        
        c = exp(hyper(M+1)); lam = hyper(M+2); rho = hyper(M+3);
        poly = hyper(1:M);
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);

        t=sum(K)/(n-1);
        
        if issym == 1
            xi = Psi2*t';
            %             XI = repmat(xi,1,N2);
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
            %             XI1 = repmat(xi1,1,N2);
            %             XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K*Psi2';
        end
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        
        
        for i = 1:M
            for j = i+1:M
                
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    TMP = ((bt1*cc1)*((bt2*cc2).*xi.^(j-i))').*(Otmp.^i);
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    O = O + c^(i+j)*(((bt1*cc1)*((bt2*cc2).*xi2.^(j-i))').*(Otmp.^i)+...
                        (((bt1j*cc2).*xi1.^(j-i))*(bt2i*cc1)').*(Otmp.^i));
                end
            end
        end
        
        if nargout > 2
            
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*K*(Psi2.*repmat((bt2*cc2).*xi.^(j-1),1,n-1))';
            end
            
        end
        
    case 'DC_eig-'
        % hp = [a1 a2 ... aM c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi Kdc Psi') + c2 (Psi Kdc Psi').^2 + ... +
        % cM (Psi Kdc Psi').^M + off-diagonal blocks
        
        c = -exp(hyper(M+1)); lam = hyper(M+2); rho = hyper(M+3);
        poly = hyper(1:M);
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        
        t=sum(K)/(n-1);
        
        if issym == 1
            xi = Psi2*t';
            %             XI = repmat(xi,1,N2);
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
            %             XI1 = repmat(xi1,1,N2);
            %             XI2 = repmat(xi2,1,N1);
            Otmp = Psi1*K*Psi2';
        end
        
        for m = 1:M
            bt1 = B1(:,1:M-m+1);
            bt2 = B2(:,1:M-m+1);
            
            cc = zeros(M-m+1,1); cc(1) = 1;
            for ii = m+1:M
                cc(ii-m+1) = cc(ii-m)*ii/(ii-m);
            end
            cc = cc.*poly(m:M);
            
            O = O + ((bt1*cc)*(bt2*cc)').*(c^2*Otmp).^m;
        end
        
        
        
        
        for i = 1:M
            for j = i+1:M
                
                bt1 = B1(:,1:M-i+1);
                bt2 = B2(:,1:M-j+1);
                
                cc1 = zeros(M-i+1,1); cc1(1) = 1;
                for ii = i+1:M
                    cc1(ii-i+1) = cc1(ii-i)*ii/(ii-i);
                end
                cc1 = cc1.*poly(i:M);
                
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                if issym == 1
                    TMP = ((bt1*cc1)*((bt2*cc2).*xi.^(j-i))').*(Otmp.^i);
                    O = O + c^(i+j)*(TMP+TMP');
                else
                    bt1j = B1(:,1:M-j+1);
                    bt2i = B2(:,1:M-i+1);
                    O = O + c^(i+j)*(((bt1*cc1)*((bt2*cc2).*xi2.^(j-i))').*(Otmp.^i)+...
                        (((bt1j*cc2).*xi1.^(j-i))*(bt2i*cc1)').*(Otmp.^i));
                end
            end
        end
        
        if nargout > 2
            
            O1h = c^2*K*(Psi2.*repmat(B2*((1:M)'.*poly),1,n-1))';
            for j = 2:M
                bt2 = B2(:,1:M-j+1);
                cc2 = zeros(M-j+1,1); cc2(1) = 1;
                for ii = j+1:M
                    cc2(ii-j+1) = cc2(ii-j)*ii/(ii-j);
                end
                cc2 = cc2.*poly(j:M);
                
                O1h  = O1h + c^(1+j)*K*(Psi2.*repmat((bt2*cc2).*xi.^(j-1),1,n-1))';
            end
            
        end
        
        
        
        
        
    case 'NEW_bd'
        % This kernel is only designed for M=2
        % hp = [a1 a2 c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi K Psi') + c2 (Psi K Psi').^2 + off-diagonal blocks
        
        c = exp(hyper(3)); lam = hyper(4); rho = hyper(5);
        poly = hyper(1:2);
        
        [TI,TJ] = meshgrid(2:n, 2:n);
        indxi = zeros(n*(n-1)/2,1);
        indxj = zeros(n*(n-1)/2,1);
        for i = 1:n-1
            ii = (2*n-1-i)*i/2;
            indxi(ii-n+1+i:ii) = diag(TI,i-1);
            indxj(ii-n+1+i:ii) = diag(TJ,i-1);
        end
        
        indx = [indxj,indxi];
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        P1 = c^2*K;
        O1 = Psi1*P1*Psi2';
        
        P2 = zeros(n*(n-1)/2,n*(n-1)/2);
        for i = 1:n*(n-1)/2
            for j = 1:n*(n-1)/2
                t = indx(i,:);
                s = indx(j,:);
                P2(i,j) = lam^(t(1)+t(2)+s(1)+s(2))*rho^(abs(t(1)-s(1)+t(2)-s(2)));
                %0.5*(lam^(t(1)+t(2)+s(1)+s(2))*rho^(abs(t(1)-s(1))+abs(t(2)-s(2)))+...
                %lam^(t(1)+t(2)+s(1)+s(2))*rho^(abs(t(1)-s(2))+abs(t(2)-s(1))));
                
            end
        end
        P2 = c^4*P2;
        
        Phi21 = zeros(N2,n*(n-1)/2);
        Phi22 = zeros(N2,n*(n-1)/2);
        for i = 1:n-1
            ii = (2*n-1-i)*i/2;
            Phi21(:,ii-n+i+1:ii) = Psi2(:,1:n-i);
            if i == 1
                Phi22(:,ii-n+i+1:ii) = Psi2;
            else
                Phi22(:,ii-n+i+1:ii) = 2*Psi2(:,i:n-1);
            end
        end
        Phi2 = Phi21.*Phi22;
        
        b12 = B2(:,1:2)*(poly.*(1:2)');
        b22 = B2(:,1)*poly(2);
        if issym == 1
            O2 = Phi2*P2*Phi2';
            O = (b12*b12').*O1 + (b22*b22').*O2;
        else
            Phi11 = zeros(N1,n*(n-1)/2);
            Phi12 = zeros(N1,n*(n-1)/2);
            for i = 1:n-1
                ii = (2*n-1-i)*i/2;
                Phi11(:,ii-n+i+1:ii) = Psi1(:,1:n-i);
                if i == 1
                    Phi12(:,ii-n+i+1:ii) = Psi1;
                else
                    Phi12(:,ii-n+i+1:ii) = 2*Psi1(:,i:n-1);
                end
            end
            Phi1 = Phi11.*Phi12;
            O2 = Phi1*P2*Phi2';
            
            b11 = B1(:,1:2)*(poly.*(1:2)');
            b21 = B1(:,1)*poly(2);
            O = (b11*b12').*O1 + (b21*b22').*O2;
        end
        
        if nargout > 2
            O1h = P1*(Psi2.*repmat(b12,1,n-1))';
        end
        
    case 'NEW_new+'
        % This kernel is only designed for M=2
        % hp = [a1 a2 c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi K Psi') + c2 (Psi K Psi').^2 + off-diagonal blocks
        
        c = exp(hyper(3)); lam = hyper(4); rho = hyper(5);
        poly = hyper(1:2);
        
        [TI,TJ] = meshgrid(2:n, 2:n);
        indxi = zeros(n*(n-1)/2,1);
        indxj = zeros(n*(n-1)/2,1);
        for i = 1:n-1
            ii = (2*n-1-i)*i/2;
            indxi(ii-n+1+i:ii) = diag(TI,i-1);
            indxj(ii-n+1+i:ii) = diag(TJ,i-1);
        end
        
        indx = [indxj,indxi];
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        P1 = c^2*K;
        O1 = Psi1*P1*Psi2';
        
        P2 = zeros(n*(n-1)/2,n*(n-1)/2);
        for i = 1:n*(n-1)/2
            for j = 1:n*(n-1)/2
                t = indx(i,:);
                s = indx(j,:);
                P2(i,j) = lam^(t(1)+t(2)+s(1)+s(2))*rho^(abs(t(1)-s(1)+t(2)-s(2)));
            end
        end
        P2 = c^4*P2;
        
        P12 = zeros(n-1,n*(n-1)/2);
        for i = 1:n-1
            for j = 1:n*(n-1)/2
                t = i+1;
                s = indx(j,:);
                P12(i,j) = lam^(t+s(1)+s(2))*rho^(abs(t-s(1)-s(2)));
            end
        end
        P12 = c^3*P12;
        
        Phi21 = zeros(N2,n*(n-1)/2);
        Phi22 = zeros(N2,n*(n-1)/2);
        for i = 1:n-1
            ii = (2*n-1-i)*i/2;
            Phi21(:,ii-n+i+1:ii) = Psi2(:,1:n-i);
            if i == 1
                Phi22(:,ii-n+i+1:ii) = Psi2;
            else
                Phi22(:,ii-n+i+1:ii) = 2*Psi2(:,i:n-1);
            end
        end
        Phi2 = Phi21.*Phi22;
        
        b12 = B2(:,1:2)*(poly.*(1:2)');
        b22 = B2(:,1)*poly(2);
        if issym == 1
            O2 = Phi2*P2*Phi2'; O12 = Psi2*P12*Phi2';
            O = (b12*b12').*O1 + (b22*b22').*O2 + (b12*b22').*O12 + (b22*b12').*O12';
        else
            Phi11 = zeros(N1,n*(n-1)/2);
            Phi12 = zeros(N1,n*(n-1)/2);
            for i = 1:n-1
                ii = (2*n-1-i)*i/2;
                Phi11(:,ii-n+i+1:ii) = Psi1(:,1:n-i);
                if i == 1
                    Phi12(:,ii-n+i+1:ii) = Psi1;
                else
                    Phi12(:,ii-n+i+1:ii) = 2*Psi1(:,i:n-1);
                end
            end
            Phi1 = Phi11.*Phi12;
            O2 = Phi1*P2*Phi2';
            
            b11 = B1(:,1:2)*(poly.*(1:2)');
            b21 = B1(:,1)*poly(2);
            
            O12 = Psi1*P12*Phi2'; O21 = Phi1*P12'*Psi2';
            O = (b11*b12').*O1 + (b21*b22').*O2 + (b11*b22').*O12 + (b21*b12').*O21;
        end
        
        if nargout > 2
            O1h = P1*(Psi2.*repmat(b12,1,n-1))' + P12*(Phi2.*repmat(b22,1,n*(n-1)/2))';
        end
        
        
    case 'NEW_new-'
        % This kernel is only designed for M=2
        % hp = [a1 a2 c lam rho sig^2 a0]
        % DC_wiener_full: c1 (Psi K Psi') + c2 (Psi K Psi').^2 + off-diagonal blocks
        
        c = -exp(hyper(3)); lam = hyper(4); rho = hyper(5);
        poly = hyper(1:2);
        
        [TI,TJ] = meshgrid(2:n, 2:n);
        indxi = zeros(n*(n-1)/2,1);
        indxj = zeros(n*(n-1)/2,1);
        for i = 1:n-1
            ii = (2*n-1-i)*i/2;
            indxi(ii-n+1+i:ii) = diag(TI,i-1);
            indxj(ii-n+1+i:ii) = diag(TJ,i-1);
        end
        
        indx = [indxj,indxi];
        
        ts = 2:n;
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        P1 = c^2*K;
        O1 = Psi1*P1*Psi2';
        
        P2 = zeros(n*(n-1)/2,n*(n-1)/2);
        for i = 1:n*(n-1)/2
            for j = 1:n*(n-1)/2
                t = indx(i,:);
                s = indx(j,:);
                P2(i,j) = lam^(t(1)+t(2)+s(1)+s(2))*rho^(abs(t(1)-s(1)+t(2)-s(2)));
            end
        end
        P2 = c^4*P2;
        
        P12 = zeros(n-1,n*(n-1)/2);
        for i = 1:n-1
            for j = 1:n*(n-1)/2
                t = i+1;
                s = indx(j,:);
                P12(i,j) = lam^(t+s(1)+s(2))*rho^(abs(t-s(1)-s(2)));
            end
        end
        P12 = c^3*P12;
        
        Phi21 = zeros(N2,n*(n-1)/2);
        Phi22 = zeros(N2,n*(n-1)/2);
        for i = 1:n-1
            ii = (2*n-1-i)*i/2;
            Phi21(:,ii-n+i+1:ii) = Psi2(:,1:n-i);
            if i == 1
                Phi22(:,ii-n+i+1:ii) = Psi2;
            else
                Phi22(:,ii-n+i+1:ii) = 2*Psi2(:,i:n-1);
            end
        end
        Phi2 = Phi21.*Phi22;
        
        b12 = B2(:,1:2)*(poly.*(1:2)');
        b22 = B2(:,1)*poly(2);
        if issym == 1
            O2 = Phi2*P2*Phi2'; O12 = Psi2*P12*Phi2';
            O = (b12*b12').*O1 + (b22*b22').*O2 + (b12*b22').*O12 + (b22*b12').*O12';
        else
            Phi11 = zeros(N1,n*(n-1)/2);
            Phi12 = zeros(N1,n*(n-1)/2);
            for i = 1:n-1
                ii = (2*n-1-i)*i/2;
                Phi11(:,ii-n+i+1:ii) = Psi1(:,1:n-i);
                if i == 1
                    Phi12(:,ii-n+i+1:ii) = Psi1;
                else
                    Phi12(:,ii-n+i+1:ii) = 2*Psi1(:,i:n-1);
                end
            end
            Phi1 = Phi11.*Phi12;
            O2 = Phi1*P2*Phi2';
            
            b11 = B1(:,1:2)*(poly.*(1:2)');
            b21 = B1(:,1)*poly(2);
            
            O12 = Psi1*P12*Phi2'; O21 = Phi1*P12'*Psi2';
            O = (b11*b12').*O1 + (b21*b22').*O2 + (b11*b22').*O12 + (b21*b12').*O21;
        end
        
        if nargout > 2
            O1h = P1*(Psi2.*repmat(b12,1,n-1))' + P12*(Phi2.*repmat(b22,1,n*(n-1)/2))';
        end
end








end


