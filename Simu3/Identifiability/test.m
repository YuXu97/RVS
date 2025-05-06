clear; n =5; N = 20; lam = rand; rho = rand; c = rand;
u = randn(N,1);
[TI,TJ] = meshgrid(2:n, 2:n);
indxi = zeros(n*(n-1)/2,1);
indxj = zeros(n*(n-1)/2,1);
for i = 1:n-1
    ii = (2*n-1-i)*i/2;
    indxi(ii-n+1+i:ii) = diag(TI,i-1);
    indxj(ii-n+1+i:ii) = diag(TJ,i-1);
end

indx = [indxj,indxi];


Psi2 = CalculatePsi(u, n);
Psi2 = Psi2(:,2:end);
Phi21 = zeros(N-n+1,n*(n-1)/2);
Phi22 = zeros(N-n+1,n*(n-1)/2);
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


P2 = zeros(n*(n-1)/2,n*(n-1)/2);
for i = 1:n*(n-1)/2
    for j = 1:n*(n-1)/2
        t = indx(i,:);
        s = indx(j,:);
        P2(i,j) = lam^(t(1)+t(2)+s(1)+s(2))*rho^(abs(t(1)-s(1)+t(2)-s(2)));
    end
end
P2 = c^4*P2;

O2 = Phi2*P2*Phi2';

% ts = 2:n;
% U = (lam*rho).^(ts)';
% V = (lam/rho).^(ts)';
% K = tril(U*V')+triu(V*U',1);
% K = c^2*K;
% 
% 
% Oo = (Psi2*K*Psi2').^2;

O3 = zeros(N-n+1,N-n+1);
for i = n:N
    for j = n:N
        for sig1 = 2:n
            for sig2 = 2:n
                for tau1 = 2:n
                    for tau2 = 2:n
                        O3(i-n+1,j-n+1) = O3(i-n+1,j-n+1) + c^4*lam^(sig1+sig2+tau1+tau2)*rho^(abs(tau1-sig1+tau2-sig2))...
                            *u(i-sig1+1)*u(i-sig2+1)*u(j-tau1+1)*u(j-tau2+1);
                    end
                end
            end
        end
    end
end

% M = 1000; n = 50;
% t = (0:M-1)';
% u = 2*rand(M,1)-1;
%
% % numerator = [2*rand-1 0];
% % denominator = [1 1.9*rand-0.95]; %[1 0.9*rand];
% % poles = -1.9*rand+0.95; %0.9*rand;
% % system = tf(numerator,denominator,1);
%
% % system = generate_linear_system_randomly(1, 0.9, 0.1);
%
% x1 = -0.9*rand;
% x2 = -0.9*rand;
% numerator = [2*rand-1 0 0];
% denominator = [1 -(x1+x2) x1*x2];
% poles = [x1 x2];
% system = tf(numerator,denominator,1);
%
% yL = lsim(system, u, t);
% [g,~] = impulse(system, t);
%
% Psi = CalculatePsi(u, n);
%
% ye = Psi*g(1:n);
% yL = yL(n:M);
%
% norm(ye - yL)
function [Psi] = CalculatePsi(x, n)
%Input: x --> input, n --> memory length
% N= length(x);
% Psi = zeros(N-n,n);
% for j = 1:n
%     Psi(:,j) = x(n-j+1:N-j);
% end
N= length(x);
Psi = zeros(N-n+1,n);

for j=1:n
    Psi(:,j) = x(n-j+1:N-j+1);
end

end
%
% function [system, poles] = generate_linear_system_randomly(order, poleub, polelb, fs_index)
% % Randomly generate a stable lienar system
% % the argument 'order' is the needed order of the system
% if nargin < 4 || strcmp(fs_index,'fast')
%     pole_max = 2;
%     pole_min = 0;
%     % pole check
%     value = 1;
% %     poles_real_min = -1;
% %     poles_tan_max = 10000;
% %
% %
%     while pole_max > poleub || pole_min < polelb  || value > 1e6
%         bw = inf;
%         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
%             mc = rss(order,1,1);
%             bw = bandwidth(mc);
%         end
%         f = 1; %bw*3*2*pi;
%         md = c2d(mc,1/f,'zoh');
%         md.d = 1;
%         system = idpoly(md);
%
%         poles = pole(system);
%         pole_max = max(abs(poles));
%         pole_min = min(abs(poles));
%
%         value = max(abs([system.f system.b]));
%     end
%
%
% %    num_con = 0;
% %
% %     while pole_max > poleub || pole_min < polelb  || value > 1e6 || num_con < ceil(order*4/5)
% %         bw = inf;
% %         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
% %             mc = rss(order,1,1);
% %             bw = bandwidth(mc);
% %         end
% %         f = 1;%bw*3*2*pi;
% %         md = c2d(mc,1/f,'zoh');
% %         md.d = 2*rand-1;
% %         system = idpoly(md);
% %
% %         poles = pole(system);
% %         pole_max = max(abs(poles));
% %         pole_min = min(abs(poles));
% %
% %         con1 = find(real(poles)<0);
% %         con2 = con1;
% %         num_con = length(intersect(con1,con2));
% %
% %         value = max(abs([system.f system.b]));
% %     end
%
%
% %    num_con = 0;
% %
% %     while pole_max > poleub || pole_min < polelb  || value > 1e6 || num_con < ceil(order)
% %         bw = inf;
% %         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
% %             mc = rss(order,1,1);
% %             bw = bandwidth(mc);
% %         end
% %         f = 1;%bw*3*2*pi;
% %         md = c2d(mc,1/f,'zoh');
% %         md.d = 2*rand-1;
% %         system = idpoly(md);
% %
% %         poles = pole(system);
% %         pole_max = max(abs(poles));
% %         pole_min = min(abs(poles));
% %
% %         pole_real = real(poles);
% %         pole_imag = imag(poles);
% %         con1 = find(pole_real>0);
% %         con2 = find(pole_imag==0)
% %         num_con = length(intersect(con1,con2));
% %
% %         value = max(abs([system.f system.b]));
% %     end
%
%
%     poles = pole(system);
%     num_poles = length(pole(system));
%     num_zeros = length(zero(system));
%     pz = [num_poles num_zeros];
%
%
% elseif strcmp(fs_index,'slow')
%
%     pole_max = 0.5; % pole check
%     value = 1;
%     while pole_max < poleub || pole_max > 1 - eps || value > 1e6
%         bw = inf;
%         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
%             mc = rss(order,1,1);
%             bw = bandwidth(mc);
%         end
%         f = bw*3*2*pi;
%         md = c2d(mc,1/f,'zoh');
%         md.d = 0;
%         system = idpoly(md);
%         pole_max = max(abs(pole(system)));
%         value = max(abs([system.f system.b]));
%     end
% end
%
% end
