function [O, U1,U2] = CalculateOutputKernel(Psi1, Psi2, M, kernel, hyper, issym)
%issym: use chol or svd to enforce the Output kernel to be symmetric (1) or
%not (0).
% in this case, Psi1 = Psi2
%O: N1 x N2 matrix
%if O is a rank-1 (or low rank) kernel matrix: O=V1*V2'
[N1, ~] = size(Psi1);
[N2, n] = size(Psi2);

O = zeros(N1,N2);
switch kernel
    case 'AMLS2os-bd'
        % hp = [c1 c2 ... cM lam rho al]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        rho = hyper(M+2);
        al = hyper(M+3);
        
        
        ts = (1:n)';
        rl = lam.^ts;
        K = (rl*rl').*cos(al*abs(ts'-ts));
        
        Otmp = Psi1*K*Psi2';
        for i = 1:M
            O = O + c(i)*Otmp.^i;
        end
        
    case 'AMLS2od-bd'
        % hp = [c1 c2 ... cM lam rho omg epsi]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        rho = hyper(M+2);
        omg = hyper(M+3);
        epsi = hyper(M+4);
        
        ts = (1:n)';
        rl = lam.^ts.*cos(omg*ts)+1+epsi;
        K = (rl*rl').*(rho.^abs(ts'-ts));
        
        Otmp = Psi1*K*Psi2';
        for i = 1:M
            O = O + c(i)*Otmp.^i;
        end
        
    case 'DC-bd'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        ts = 1:n;
        c = exp(hyper(1:M));
        lam = hyper(M+1);
        rho = hyper(M+2);
        rl = lam.^ts';
        K = (rl*rl').*(rho.^abs(ts'-ts));
        
        %         if issym == 1
        %             try
        %                 Kc = chol(K)';
        %                 tmp = Psi2*Kc;
        %             catch
        %                 [Uk,Sk,~] = svd(K);
        %                 tmp = Psi2*Uk*sqrt(Sk);
        %             end
        %             Otmp = tmp*tmp';
        %         else
        %             Otmp = Psi1*K*Psi2';
        %         end
        Otmp = Psi1*K*Psi2';
        for i = 1:M
            O = O + c(i)*Otmp.^i;
        end
        

    case 'DC-bd-odd'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        ts = 1:n;
        pl = hyper(1:M);
        c = exp(hyper(M+1));
        lam = hyper(M+2);
        rho = hyper(M+3);
        rl = lam.^ts';
        K = (rl*rl').*(rho.^abs(ts'-ts));
                %         if issym == 1
        %             try
        %                 Kc = chol(K)';
        %                 tmp = Psi2*Kc;
        %             catch
        %                 [Uk,Sk,~] = svd(K);
        %                 tmp = Psi2*Uk*sqrt(Sk);
        %             end
        %             Otmp = tmp*tmp';
        %         else
        %             Otmp = Psi1*K*Psi2';
        %         end
        
        Otmp = c*Psi1*K*Psi2';
%         G = load(['g.mat']);
%         g = G.g;
%         Otmp = c*(Psi1*g)*(Psi2*g)';
        
%         for i = 1:M
%             O = O + pl(i)^2*Otmp.^(2*i-1);
%         end
        
%         O = zeros(N1,N2);
        TMP = Otmp;
        DP = Otmp.^2;
        for i = 1:M
            O = O + pl(i)^2*TMP;
            TMP = TMP.*DP;
        end       
    
        
    case '2DC_SI2od'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        ts = (1:n)';
        c1 = exp(hyper(1));
        c2 = exp(hyper(2));
        lam1 = hyper(3);
        lam2 = hyper(4);
        rho1 = hyper(5);
        rho2 = hyper(6);
        
        r1 = lam1.^ts;
        r2 = lam2.^ts;
        K1 = (r1*r1').*(rho1.^abs(ts'-ts));
        K2 = (r2*r2').*(rho2.^abs(ts'-ts));
        
        lama = hyper(7);
        xi = hyper(8);
        lam = hyper(9);
        
        al = -log(lama);
        bt = al/xi*sqrt(1-xi^2);
        gm = -log(lam);
        amg = al-gm;
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
        
        
        O = Psi1*(c1*K1+c2*K2+K)*Psi2';
        
        
     case '2DC_SI2od-bd-odd-polyfix'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        pl = [-0.898538200392490,1.53751835545513,-0.771214612369874,0.219960559272809,-0.0372181216502121,0.00380527430788070,-0.000230461565158972,7.60227506719557e-06,-1.05164197121184e-07];
        ts = (1:n)';
        c1 = exp(hyper(1));
        c2 = exp(hyper(2));
        lam1 = hyper(3);
        lam2 = hyper(4);
        rho1 = hyper(5);
        rho2 = hyper(6);
        
        r1 = lam1.^ts;
        r2 = lam2.^ts;
        K1 = (r1*r1').*(rho1.^abs(ts'-ts));
        K2 = (r2*r2').*(rho2.^abs(ts'-ts));
        
        lama = hyper(7);
        xi = hyper(8);
        lam = hyper(9);
        
        al = -log(lama);
        bt = al/xi*sqrt(1-xi^2);
        gm = -log(lam);
        amg = al-gm;
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
        
        
        Otmp = Psi1*(c1*K1+c2*K2+K)*Psi2';
        
        TMP = Otmp;
        DP = Otmp.^2;
        for i = 1:M
            O = O + pl(i)^2*TMP;
            TMP = TMP.*DP;
        end
        
        
     case 'SI2od_dc-bd'
        ts = (1:n)';
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
        
        for i = 1:M
            O = O + pl(i)^2*Otmp.^i;
        end
        
        
     case 'SI2od_dc-bd-odd-polyfix'
        %pl = [-0.898538200392490,1.53751835545513,-0.771214612369874,0.219960559272809,-0.0372181216502121,0.00380527430788070,-0.000230461565158972,7.60227506719557e-06,-1.05164197121184e-07];
        pl = [-0.987735385675577,1.91515217210649,-1.24354495221076,0.491490093969614,-0.123171525837517,0.0201820981177221,-0.00218797053322340,0.000155425229728045,-6.95085209452765e-06,1.77418233001695e-07,-1.96993920710876e-09];
        ts = (1:n)';
        c = exp(hyper(1));

        lam = hyper(2);
        rho = hyper(3);
        
        lama = hyper(4);
        xi = hyper(5);
        lamg = hyper(6);
        
        al = -log(lama);
        bt = al/xi*sqrt(1-xi^2);
        gm = -log(lamg);
        amg = al-gm;
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
        
        TMP = Otmp;
        DP = Otmp.^2;
        for i = 1:M
            O = O + pl(i)^2*TMP;
            TMP = TMP.*DP;
        end
        
    case 'SI2od'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        ts = (1:n)';
        c = exp(hyper(1));

        lama = hyper(2);
        xi = hyper(3);
        lam = hyper(4);
        
        al = -log(lama);
        bt = al/xi*sqrt(1-xi^2);
        gm = -log(lam);
        amg = al-gm;
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
        
        
        O = c*Psi1*K*Psi2';
        
    case 'SI2od_dc'
        ts = (1:n)';
        c = exp(hyper(1));
        
        lam = hyper(2);
        rho = hyper(3);
        
        lama = hyper(4);
        xi = hyper(5);
        lamg = hyper(6);
        
        al = -log(lama);
        bt = al/xi*sqrt(1-xi^2);
        gm = -log(lamg);
        amg = al-gm;
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
           
        O = Psi1*(c*Kdc+K)*Psi2';
        
    case 'SI2od_dc-bd-odd'
        pl = hyper(1:M);
        ts = (1:n)';
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
        
        TMP = Otmp;
        DP = Otmp.^2;
        for i = 1:M
            O = O + pl(i)^2*TMP;
            TMP = TMP.*DP;
        end
        
        
     
    
    case 'SI2od_tdc'
        ts = (0:n-1)';
        c = exp(hyper(1));
        
        lam = hyper(2);
        rho = hyper(3);
        
        lama = hyper(4);
        xi = hyper(5);
        lamg = hyper(6);
        
        al = -log(lama);
        bt = al/xi*sqrt(1-xi^2);
        gm = -log(lamg);
        amg = al-gm;
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
        
        rld = ts.*lam.^ts;
        Kdc = (rld*rld').*(rho.^abs(ts'-ts));
           
        O = Psi1*(c*Kdc+K)*Psi2';
        
        
    case 'hDC-bd-odd-polyfix'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        pl = [-0.987735385675577,1.91515217210649,-1.24354495221076,0.491490093969614,-0.123171525837517,0.0201820981177221,-0.00218797053322340,0.000155425229728045,-6.95085209452765e-06,1.77418233001695e-07,-1.96993920710876e-09];
       
        ts = 1:n;

        c1 = exp(hyper(1));
        lam = hyper(2);
        rho = hyper(3);
        lamn = hyper(4);
        var1 = hyper(5);
        var2 = hyper(6);
        omg = hyper(7);
        rl = lam.^ts';
        rl2 = lamn.^ts';
        gm1 = (var1+var2)/2;
        gm2 = (var1-var2)/2;
        K = c1*(rl*rl').*(rho.^abs(ts'-ts))+(rl2*rl2')...
            .*(gm1*cos(omg*pi*(ts'-ts))+gm2*cos(omg*pi*(ts'+ts)));
        Otmp = Psi1*K*Psi2';
        

        TMP = Otmp;
        DP = Otmp.^2;
        for i = 1:M
            O = O + pl(i)^2*TMP;
            TMP = TMP.*DP;
        end
        
    case 'hDC_2dc-bd-odd-polyfix'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        pl = [-0.987735385675577,1.91515217210649,-1.24354495221076,0.491490093969614,-0.123171525837517,0.0201820981177221,-0.00218797053322340,0.000155425229728045,-6.95085209452765e-06,1.77418233001695e-07,-1.96993920710876e-09];
        
        ts = 1:n;

        c1 = exp(hyper(1));
        c2 = exp(hyper(2));
        lam1 = hyper(3);
        rho1 = hyper(4);
        lam2 = hyper(5);
        rho2 = hyper(6);
        
        lamn = hyper(7);
        var1 = hyper(8);
        var2 = hyper(9);
        omg = hyper(10);
        r1 = lam1.^ts';
        r2 = lam2.^ts';
        
        rl2 = lamn.^ts';
        gm1 = (var1+var2)/2;
        gm2 = (var1-var2)/2;
        K = c1*(r1*r1').*(rho1.^abs(ts'-ts))+c2*(r2*r2').*(rho2.^abs(ts'-ts))...
            +(rl2*rl2').*(gm1*cos(omg*pi*(ts'-ts))+gm2*cos(omg*pi*(ts'+ts)));
        Otmp = Psi1*K*Psi2';
        

        TMP = Otmp;
        DP = Otmp.^2;
        for i = 1:M
            O = O + pl(i)^2*TMP;
            TMP = TMP.*DP;
        end
    
    case 'hDC-bd'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        ts = 0:n-1;

        pl = hyper(1:M);
        c1 = exp(hyper(M+1));
        lam = hyper(M+2);
        rho = hyper(M+3);
        lamn = hyper(M+4);
        var1 = hyper(M+5);
        var2 = hyper(M+6);
        omg = hyper(M+7);
        rl = lam.^ts';
        rl2 = lamn.^ts';
        gm1 = (var1+var2)/2;
        gm2 = (var1-var2)/2;
        K = c1*(rl*rl').*(rho.^abs(ts'-ts))+(rl2*rl2')...
            .*(gm1*cos(omg*pi*(ts'-ts))+gm2*cos(omg*pi*(ts'+ts)));
        

        Otmp = Psi1*K*Psi2';
        for i = 1:M
            O = O + pl(i)^2*Otmp.^i;
        end
        
    case 'hDC-bd-odd'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        ts = 0:n-1;

        pl = hyper(1:M);
        c1 = exp(hyper(M+1));
        lam = hyper(M+2);
        rho = hyper(M+3);
        lamn = hyper(M+4);
        var1 = hyper(M+5);
        var2 = hyper(M+6);
        omg = hyper(M+7);
        rl = lam.^ts';
        rl2 = lamn.^ts';
        gm1 = (var1+var2)/2;
        gm2 = (var1-var2)/2;
        K = c1*(rl*rl').*(rho.^abs(ts'-ts))+(rl2*rl2')...
            .*(gm1*cos(omg*pi*(ts'-ts))+gm2*cos(omg*pi*(ts'+ts)));
        Otmp = Psi1*K*Psi2';
        

        TMP = Otmp;
        DP = Otmp.^2;
        for i = 1:M
            O = O + pl(i)^2*TMP;
            TMP = TMP.*DP;
        end
        
     case 'hOS'

        ts = 0:n-1;
        c1 = exp(hyper(1));
        lamn = hyper(2);
        var1 = hyper(3);
        var2 = hyper(4);
        omg = hyper(5);
        rl2 = lamn.^ts';
        gm1 = (var1+var2)/2;
        gm2 = (var1-var2)/2;
        K = c1*(rl2*rl2').*(gm1*cos(omg*pi*(ts'-ts))+gm2*cos(omg*pi*(ts'+ts)));
        O = Psi1*K*Psi2';
        
% 
%         TMP = Otmp;
%         DP = Otmp.^2;
%         for i = 1:M
%             O = O + pl(i)^2*TMP;
%             TMP = TMP.*DP;
%         end
        
    case 'DC-bd-odd-polyfix'
        % hp = [c1 c2 ... cM lam rho]
        % DC-bd: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M
        % cm = am^2 c^(2m)
        ts = 0:n-1;
        %c = exp(hyper(1));
        %lam = hyper(2);
        %rho = hyper(3);
        %rl = ts'.*lam.^ts';
        %K = (rl*rl').*(rho.^abs(ts'-ts));
        %Otmp = c*Psi1*K*Psi2';
        
        c1 = exp(hyper(1));
        lam = hyper(2);
        rho = hyper(3);
        lamn = hyper(4);
        var1 = hyper(5);
        var2 = hyper(6);
        omg = hyper(7);
        rl = lam.^ts';
        rl2 = lamn.^ts';
        gm1 = (var1+var2)/2;
        gm2 = (var1-var2)/2;
        K = c1*(rl*rl').*(rho.^abs(ts'-ts))+(rl2*rl2')...
            .*(gm1*cos(omg*pi*(ts'-ts))+gm2*cos(omg*pi*(ts'+ts)));
        Otmp = Psi1*K*Psi2';
        
        
        pl = [2.097474850579006,-1.510904700754373,0.524293124181930,-0.062441127759624];
        TMP = Otmp;
        DP = Otmp.^2;
        for i = 1:M
            O = O + pl(i)^2*TMP;
            TMP = TMP.*DP;
        end
        
    case 'DC-bd-odd-r1'
        g = load('g.mat').g;
        pl = [2.34304190766050,-2.51441753277124,1.67027195631447,...
            -0.591409073750925,0.104839744000900,-0.00732245172673269];          
        eta1 = Psi1*g;
        eta2 = Psi2*g;
        
        U1 = pl.*eta1.^(1:2:(2*M-1));
        U2 = pl.*eta2.^(1:2:(2*M-1));
       
        
        
        
    case 'DC-opt'
        % hp = [c1 c2 ... cM lam rho]
        % DC-opt: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
        % cm = am c^m
        
        ts = 1:n;
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        rho = hyper(M+2);
        rl = lam.^ts';
        K = (rl*rl').*(rho.^abs(ts'-ts));
        t = (lam*rho).^(ts);
        
        if issym == 1
            xi = Psi2*t';
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
            Otmp = Psi1*K*Psi2';
        end
        
        for i = 1:M
            O = O + c(i)^2*Otmp.^i;
        end
        
        
        for i = 1:M
            for j = i+1:M
                if issym == 1
                    tmpp = (xi.^i)*(xi.^j)';
                    O = O + c(i)*c(j)*(tmpp+tmpp');
                else
                    O = O + c(i)*c(j)*((xi1.^i)*(xi2.^j)'+(xi1.^j)*(xi2.^i)');
                end
            end
        end
        
        
        
        
    case 'DC-dc'
        % hp = [c1 c2 ... cM lam rho]
        % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
        % cm = am c^m
        
        ts = 1:n;
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        rho = hyper(M+2);
        rl = lam.^ts';
        K = (rl*rl').*(rho.^abs(ts'-ts));
        
        
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
        
    case 'DC-dcp'
        % hp = [c1 c2 ... cM lam rho]
        % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
        % cm = am c^m
        
        ts = 1:n;
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        rho = hyper(M+2);
        
        
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        t = (lam*rho).^(ts)*sqrt(1-rho^2)/rho;
        
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
                if issym == 1
                    TMP = (Otmp.^i).*(XI.^(j-i));
                    O = O + c(i)*c(j)*(TMP+TMP');
                else
                    O = O + c(i)*c(j)*((Otmp.^i).*(XI2.^(j-i))'+(Otmp.^i).*(XI1.^(j-i)));
                end
            end
        end
        
        
    case 'DC-dcr'
        % hp = [c1 c2 ... cM lam rho alpha]
        % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
        % cm = am c^m
        
        ts = 1:n;
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        rho = hyper(M+2);
        alpha = hyper(M+3);
        
        
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        t = alpha*(lam*rho).^(ts)*sqrt(1-rho^2)/rho;
        
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
                if issym == 1
                    TMP = (Otmp.^i).*(XI.^(j-i));
                    O = O + c(i)*c(j)*(TMP+TMP');
                else
                    O = O + c(i)*c(j)*((Otmp.^i).*(XI2.^(j-i))'+(Otmp.^i).*(XI1.^(j-i)));
                end
            end
        end
        
    case 'DC-dcs'
        % hp = [c1 c2 ... cM lam rho alpha]
        % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
        % cm = am c^m
        
        ts = 1:n;
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        rho = hyper(M+2);
        alpha = hyper(M+3);
        
        
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        t = (alpha*sqrt(1-rho^2)/rho+1-alpha)*(lam*rho).^(ts);
        
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
                if issym == 1
                    TMP = (Otmp.^i).*(XI.^(j-i));
                    O = O + c(i)*c(j)*(TMP+TMP');
                else
                    O = O + c(i)*c(j)*((Otmp.^i).*(XI2.^(j-i))'+(Otmp.^i).*(XI1.^(j-i)));
                end
            end
        end
        
        
    case 'DC-ob'
        % hp = [c1 c2 ... cM lam rho]
        
        
        ts = 1:n;
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        rho = hyper(M+2);
        rl = lam.^ts';
        K = (rl*rl').*(rho.^abs(ts'-ts));
        
        r = 100;
        d = 1./(((1:r)'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = (1:r)'-0.5;
        Ut = sqrt(2)*sin(v1*v2').*repmat((lam/rho).^ts',1,r);
        cf = sqrt(2)*sqrt(d);
        
        t = Ut*diag(sqrt(d))*cf; t = t';
        
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
        
    case 'DC-ob-ex'
        % hp = [c1 c2 ... cM lam rho alpha beta]
        
        
        ts = 1:n;
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        rho = hyper(M+2);
        %         alpha = hyper(M+3);
        %         beta = hyper(M+4);
        
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        
        r = 100;
        d = 1./(((1:r)'-0.5).^2)/pi^2;
        v1 = pi*rho.^(2*ts)'; v2 = (1:r)'-0.5;
        Ut = sqrt(2)*sin(v1*v2').*repmat((lam/rho).^ts',1,r);
        dsq = sqrt(d); signs = (-1).^[(1:r)'+1];
        cf = sqrt(2)*dsq.*signs;
        t = Ut*diag(sqrt(d))*cf; t = t';
        
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
                if issym == 1
                    TMP = (Otmp.^i).*(XI.^(j-i));
                    O = O + c(i)*c(j)*(TMP+TMP');
                else
                    O = O + c(i)*c(j)*((Otmp.^i).*(XI2.^(j-i))'+(Otmp.^i).*(XI1.^(j-i)));
                end
            end
        end
        
        
    case 'DC-eig'
        % hp = [c1 c2 ... cM lam rho alpha]
        % DC-dc: a1^2c^2 (Psi Kdc Psi') + a2^2c^4 (Psi Kdc Psi').^2 + ... + aM^2c^(2M) (Psi Kdc Psi').^M + off-diagonal blocks
        % cm = am c^m
        
        ts = 1:n;
        c = hyper(1:M).^3;
        lam = hyper(M+1);
        rho = hyper(M+2);
        
        
        U = (lam*rho).^(ts)';
        V = (lam/rho).^(ts)';
        K = tril(U*V')+triu(V*U',1);
        
        %         [Uk,Sk,~] = svd(K);
        %         d = diag(Sk); dsq = sqrt(d);
        %         t = Uk*d/norm(dsq); t = t';
        
        [V,D] = eig(K);
        d = diag(D); dsq = real(sqrt(d));
        t = V*d/norm(dsq); t = t';
        
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