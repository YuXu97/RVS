function [sys0, data] = testsystem_wienergp_4d(T, par, repi)
% -- 4th order system

T = 2*T;

sys0.nx = 4;
if(par.input), sys0.nu = 1; else sys0.nu = 0; end;
sys0.ny = 1;
sys0.type = 'Wiener';
sys0.which = '4d atan';

% Dynamics
a= [0.368, 0.888, 0.524, 0.555];
c = [1, 0.1, -0.49, 0.01];

sys = tf([0 a], [1 c], 1, 'Variable','z^-1');


sys = ss(sys);


sys0.ss.A = sys.a;
sys0.ss.C = sys.c;
if(par.input)
    sys0.ss.B = sys.b;
    u = 0.5*randn(1,T);
else
    sys0.ss.B = [];
    u = [];
end


sys0.ss.h = @(z)nonmon(z,0.5,0.3);
sys0.ss.Q = 0;%0.025^2*eye(sys0.nx);
sys0.ss.R = 0.1^2; % Measurement noise


% Initial state -----------------------------------------------------------
sys0.ss.X1 = zeros(sys0.nx,1);
sys0.ss.P1 = 1*eye(sys0.nx);

% Initialize
data.tVec = (1:T);
x = zeros(sys0.nx, T);
z = zeros(1, T);
y = zeros(1, T);

% w = mvnrnd(zeros(1,sys0.nx), sys0.ss.Q, T)';
e = sqrt(sys0.ss.R)*randn(1, T);

x(:,1) = mvnrnd(sys0.ss.X1, sys0.ss.P1, 1)';

for(t = 1:T-1)
    if(par.input)
        x(:,t+1) = sys0.ss.A*x(:,t) + sys0.ss.B*u(t);
        %x(:,t+1) = sys0.ss.A*x(:,t) + sys0.ss.B*u(t) + w(:,t);
    else
        x(:,t+1) = sys0.ss.A*x(:,t);
        %x(:,t+1) = sys0.ss.A*x(:,t) + w(:,t);
    end
    z(t) = sys0.ss.C*x(:,t);
    y(t) = sys0.ss.h(z(t)) + e(t);
    %   y(t) = sys0.ss.h(z(t));
end

z(t+1) = sys0.ss.C*x(:,t+1);
y(t+1) = sys0.ss.h(z(t+1)) + e(t+1);

N = T/2;
data.u = u(1:N);
data.uv = u(N+1:2*N);
data.x = x;
data.z = z;
data.y = y(1:N);
data.yv = y(N+1:2*N);
% data.w = w;
data.e = e;

save(['databank/' 'data_' int2str(repi) '.mat'], 'data');
save(['databank/' 'sys0_' int2str(repi) '.mat'], 'sys0');
end

