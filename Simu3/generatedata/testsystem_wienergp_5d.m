function [sys0, data] = testsystem_wienergp_5d(T,par)

load sys5d;

% -- 5th order system from Wiener-paper
sys0.nx = 5;
if(par.input), sys0.nu = 1; else sys0.nu = 0; end;
sys0.ny = 1;
sys0.type = 'Wiener';
sys0.which = '5d Wiener';

if(par.input)
    u = 3*randn(1,T);
else
    sys0.ss.B = [];
    u = [];
end

sys0.ss.h = @(z)(2*max(min(z,0.5),-0.5));     % saturation
sys0.ss.Q = 0.1^2*eye(sys0.nx);
sys0.ss.R = 0.05^2; % Measurement noise

% Initial state -----------------------------------------------------------
sys0.ss.X1 = zeros(sys0.nx,1);
sys0.ss.P1 = 1*eye(sys0.nx);

% Initialize
data.tVec = (1:T);
x = zeros(sys0.nx, T);
z = zeros(1, T);
y = zeros(1, T);

w = mvnrnd(zeros(1,sys0.nx), sys0.ss.Q, T)';
e = sqrt(sys0.ss.R)*randn(1, T);
x(:,1) = mvnrnd(sys0.ss.X1, sys0.ss.P1, 1)';

for(t = 1:T-1)
    if(par.input)
        x(:,t+1) = sys0.ss.A*x(:,t) + sys0.ss.B*u(t) + w(:,t);
    else
        x(:,t+1) = sys0.ss.A*x(:,t) + w(:,t);
    end
    z(t) = sys0.ss.C*x(:,t);
    y(t) = sys0.ss.h(z(t)) + e(t);
end
z(t+1) = sys0.ss.C*x(:,t+1);
y(t+1) = sys0.ss.h(z(t+1)) + e(t+1);


data.u = u;
data.x = x;
data.z = z;
data.y = y;
data.w = w;
data.e = e;

