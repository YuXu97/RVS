function [sys0, data] = testsystem_wienergp_2d(T, par, repi)
% -- 4th order system

T = 2*T;

sys0.nx = 2;
if(par.input), sys0.nu = 1; else sys0.nu = 0; end;
sys0.ny = 1;
sys0.type = 'Wiener';
sys0.which = '2d atan';

% Dynamics
% numerator = [0.06 0];
% denominator = [1 -1.8036 0.8338];
% sys = tf(numerator,denominator,1);
sys = generate_linear_system_randomly(2,0.95,0.1);
%
sys = ss(sys);


sys0.ss.A = sys.a;
sys0.ss.C = sys.c;
if(par.input)
    sys0.ss.B = sys.b;
    u = randn(T,1);
%     rng(2,'twister')
    Syst_filter = generate_linear_system_randomly(5,0.95,0.1);
    u = sim(Syst_filter,u);
else
    sys0.ss.B = [];
    u = [];
end
%\
cp = 2*randn(2,1)-1;
% sys0.ss.h = @(z)(2*max(min(z,0.5),-0.5));
sys0.ss.h = @(z) cp(1)*z.^2 + cp(2)*z;
sys0.ss.Q = 0;%0.025^2*eye(sys0.nx);
sys0.ss.R = 0.1^2; % Measurement noise


% Initial state -----------------------------------------------------------
sys0.ss.X1 = zeros(sys0.nx,1);
sys0.ss.P1 = 1*eye(sys0.nx);

% Initialize
data.tVec = (1:T);
x = zeros(sys0.nx, T);
warning off
z = lsim(sys, u, 1:T);

% z = zeros(1, T);
% y = zeros(1, T);

% w = mvnrnd(zeros(1,sys0.nx), sys0.ss.Q, T)';


% x(:,1) = mvnrnd(sys0.ss.X1, sys0.ss.P1, 1)';
%
% for(t = 1:T-1)
%     if(par.input)
%         x(:,t+1) = sys0.ss.A*x(:,t) + sys0.ss.B*u(t);
%         %x(:,t+1) = sys0.ss.A*x(:,t) + sys0.ss.B*u(t) + w(:,t);
%     else
%         x(:,t+1) = sys0.ss.A*x(:,t);
%         %x(:,t+1) = sys0.ss.A*x(:,t) + w(:,t);
%     end
%     z(t) = sys0.ss.C*x(:,t);
%     y(t) = sys0.ss.h(z(t));
% end
%
% z(t+1) = sys0.ss.C*x(:,t+1);
% y(t+1) = sys0.ss.h(z(t+1));

snr = db2pow(5);
ytrue = sys0.ss.h(z);
e = randn(T,1);
e = (e - mean(e))*std(ytrue)/std(e)/sqrt(snr);
y = ytrue + e;

N = T/2;
data.u = u(1:N)';
data.uv = u(N+1:2*N)';
data.x = x;
data.z = z';
data.y = y(1:N)';
data.yv = y(N+1:2*N)';
% data.w = w;
data.ytrue = ytrue;
data.e = e';

save(['databank/' 'data_' int2str(repi) '.mat'], 'data');
save(['databank/' 'sys0_' int2str(repi) '.mat'], 'sys0');
end

function [system, poles] = generate_linear_system_randomly(order, poleub, polelb, fs_index)
% Randomly generate a stable lienar system
% the argument 'order' is the needed order of the system
if nargin < 4 || strcmp(fs_index,'fast')
    pole_max = 2;
    pole_min = 0;
    % pole check
    value = 1;
    %     poles_real_min = -1;
    %     poles_tan_max = 10000;

    % ratio = 0;
    while pole_max > poleub || pole_min < polelb  || value > 1e6 %|| ratio < 0.8
        bw = inf;
        while isinf(bw)||isnan(bw) % finite bandwidth of generated system
            mc = rss(order,1,1);
            bw = bandwidth(mc);
        end
        f = 1; %bw*3*2*pi;
        md = c2d(mc,1/f,'zoh');
        md.d = 0;
        system = idpoly(md);

        poles = pole(system);
        pole_max = max(abs(poles));
        pole_min = min(abs(poles));

        value = max(abs([system.f system.b]));
        %         ratio = cos(angle(poles(1)));

    end
    %

    % %poles that are close to the real line
    % pole_max = 2;
    % pole_min = 0;
    % poles_imag_max = 2;
    %
    % while pole_max > poleub || pole_min < polelb  || value > 1e6 || poles_imag_max > 0.1
    %     bw = inf;
    %     while isinf(bw)||isnan(bw) % finite bandwidth of generated system
    %         mc = rss(order,1,1);
    %         bw = bandwidth(mc);
    %     end
    %     f = 1; %bw*3*2*pi;
    %     md = c2d(mc,1/f,'zoh');
    %     md.d = 0;
    %     system = idpoly(md);
    %
    %     poles = pole(system);
    %     pole_max = max(abs(poles));
    %     pole_min = min(abs(poles));
    %     poles_imag_max = max(imag(poles));
    %
    %     value = max(abs([system.f system.b]));
    % end
    % system = zpk([zero(system);0],pole(system),2*rand-1,1);



    % five poles with the largest modulus is real and around 0.75
    % pole_max = 100;
    %
    % while pole_max > poleub || pole_min < polelb  || value > 1e6 || pole_max > 0.45
    %     bw = inf;
    %     while isinf(bw)||isnan(bw) % finite bandwidth of generated system
    %         mc = rss(order-5,1,1);
    %         bw = bandwidth(mc);
    %     end
    %     f = 1;%bw*3*2*pi;
    %     md = c2d(mc,1/f,'zoh');
    %     md.d = 0;
    %     system_t = idpoly(md);
    %     poles = pole(system_t);
    %     pole_max = max(abs(poles));
    %     pole_min = min(abs(poles));
    %     value = max(abs([system_t.f system_t.b]));
    % end
    % zeros_zpk = [zero(system_t);rand;rand(5,1)];%make sure there is zeros delay
    % poles_zpk = [poles; 0.25*rand(5,1)+0.7];
    % system = zpk(zeros_zpk,poles_zpk,2*rand-1,1);





    %    num_con = 0;
    %
    %     while pole_max > poleub || pole_min < polelb  || value > 1e6 || num_con < ceil(order*4/5)
    %         bw = inf;
    %         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
    %             mc = rss(order,1,1);
    %             bw = bandwidth(mc);
    %         end
    %         f = 1;%bw*3*2*pi;
    %         md = c2d(mc,1/f,'zoh');
    %         md.d = 2*rand-1;
    %         system = idpoly(md);
    %
    %         poles = pole(system);
    %         pole_max = max(abs(poles));
    %         pole_min = min(abs(poles));
    %
    %         con1 = find(real(poles)<0);
    %         con2 = con1;
    %         num_con = length(intersect(con1,con2));
    %
    %         value = max(abs([system.f system.b]));
    %     end


    %    num_con = 0;
    %
    %     while pole_max > poleub || pole_min < polelb  || value > 1e6 || num_con < ceil(order)
    %         bw = inf;
    %         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
    %             mc = rss(order,1,1);
    %             bw = bandwidth(mc);
    %         end
    %         f = 1;%bw*3*2*pi;
    %         md = c2d(mc,1/f,'zoh');
    %         md.d = rand;
    %         system = idpoly(md);
    %
    %         poles = pole(system);
    %         pole_max = max(abs(poles));
    %         pole_min = min(abs(poles));
    %
    %         pole_real = real(poles);
    %         pole_imag = imag(poles);
    %         con1 = find(pole_real>0);
    %         con2 = find(pole_imag==0)
    %         num_con = length(intersect(con1,con2));
    %
    %         value = max(abs([system.f system.b]));
    %     end


    poles = pole(system);
    num_poles = length(pole(system));
    num_zeros = length(zero(system));
    pz = [num_poles num_zeros];


elseif strcmp(fs_index,'slow')

    pole_max = 0.5; % pole check
    value = 1;
    while pole_max < poleub || pole_max > 1 - eps || value > 1e6
        bw = inf;
        while isinf(bw)||isnan(bw) % finite bandwidth of generated system
            mc = rss(order,1,1);
            bw = bandwidth(mc);
        end
        f = bw*3*2*pi;
        md = c2d(mc,1/f,'zoh');
        md.d = 0;
        system = idpoly(md);
        pole_max = max(abs(pole(system)));
        value = max(abs([system.f system.b]));
    end
end

end

