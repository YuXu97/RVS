function [sys0, data] = testsystem_gen6d(T, par, repi)
% -- 4th order system

T = 4*T;

sys0.nx = 6;
if(par.input), sys0.nu = 1; else sys0.nu = 0; end;
sys0.ny = 1;
sys0.type = 'Wiener';
sys0.which = '2d atan';
% 
% % Dynamics
% Lorder = 10;
% sys = generate_linear_system_randomly(Lorder, 0.9, 0.1);


a = [-0.467, 1.12, -0.925, 0.308, -0.0364, 0.00110];
c = [-2.67, 2.96, -2.01, 0.914, -0.181, -0.0102];


sys = tf([0 a], [1 c], 1, 'Variable','z^-1'); 
g = impulse(sys,0:T-1);
sys = ss(sys);


% sys = tf([0 a]./g(1), [1 c], 1, 'Variable','z^-1');
% sys = ss(sys);
% g = impulse(sys,0:T-1);
%  
% [zs,gain] = zero(sys);
% id =find(zs>1);
% zs(id) = 1/zs(id);
% ps = pole(sys);
% sys = zpk(zs,ps,gain,1);
% 
% Dynamics
% c= [0.368, 0.888, 0.524, 0.555];
% a = [1, 0.1, -0.49, 0.01];
% 
% sys = tf([0 c], [1 a], 1, 'Variable','z^-1');

% 
% numerator = [0.06 0];
% denominator = [1 -1.8036 0.8338];
% sys = tf(numerator,denominator,1);
% sys = ss(sys);

% numerator = [0.06 0];
% denominator = [1 -0.8036];
% sys = tf(numerator,denominator,1);
% sys = ss(sys);

% numerator = [0.05 0];
% denominator = [1 -1.8206 0.8292];
% sys = tf(numerator,denominator,1);
% sys = ss(sys);

% numerator = [0.2 0];
% denominator = [1 -1.3094 0.4329];
% sys = tf(numerator,denominator,1);
% sys = ss(sys);

% numerator = [0.5 0];
% denominator = [1 -0.5338];
% sys = tf(numerator,denominator,1);
% sys = ss(sys);



sys0.ss.A = sys.a;
sys0.ss.C = sys.c;
if(par.input)
    sys0.ss.B = sys.b;
%     rng(repi,'twister')
    u = randn(T,1);
else
    sys0.ss.B = [];
    u = [];
end

% sys0.ss.h = @(z)nonmon(z,0.5,0.3);
sys0.ss.h = @(z)(2*max(min(z,0.5),-0.5));

%    sys0.ss.h = @(z) 2*z.^5-3*z.^3 + z;
% plc = 2*rand(3,1) - 1;
% 
% sys0.ss.h = @(z) plc(1)*z.^2 + plc(2)*z + plc(3);

% sys0.ss.h = @(z) 0.5*z.^2 + 0.1*z - 0.1;
% sys0.ss.h = @(z) sin(pi*z)./(pi*z);
% sys0.ss.h = @(z) z./(z+exp(-z));
sys0.ss.Q = 0;%0.025^2*eye(sys0.nx);
sys0.ss.R = 0.1^2; % Measurement noise


% Initial state -----------------------------------------------------------
sys0.ss.X1 = zeros(sys0.nx,1);
sys0.ss.P1 = 1*eye(sys0.nx);

% Initialize
data.tVec = (1:T);
x = zeros(sys0.nx, T);
warning off
z = lsim(sys, u, 0:T-1); 
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
% z(t+1) = sys0.ss.C*x(:,t+1);0;
% y(t+1) = sys0.ss.h(z(t+1));
% 
ytrue = sys0.ss.h(z);
% C = load(['coef.mat']);
% coef = C.coef;
% ytrue = zeros(T,1);
% for i = 1:6
%     ytrue = ytrue + coef(i)*z.^(2*i-1);
% end
% ytrue(1) = 1;

e = sqrt(sys0.ss.R)*randn(T,1);
%  snr = db2pow(5);
%  e = randn(T,1);
%  e = (e - mean(e))*std(ytrue)/std(e)/sqrt(snr);
y = ytrue + e;

N = T/4;

data.u = u(1:N)';
data.uv = u(N+1:2*N)';
data.ut = u(2*N+1:4*N)';
data.gtrue = g;
data.ytrue = ytrue(1:N)';
data.yvtrue = ytrue(N+1:2*N)';
data.yttrue = ytrue(2*N+1:4*N)';
data.z = z;
data.y = y(1:N)';
data.yv = y(N+1:2*N)';
data.yt = y(2*N+1:4*N)';
% data.plc = plc;
data.e = e';
data.sys = sys;
data.nonl = sys0.ss.h;

save(['databank/' 'data_' int2str(repi) '.mat'], 'data');
save(['databank/' 'sys0_' int2str(repi) '.mat'], 'sys0');
end

