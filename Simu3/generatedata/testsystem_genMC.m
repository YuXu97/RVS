function [sys0, data] = testsystem_genMC(T, par, repi)
% -- 4th order system

T = 4*T;

sys0.nx = 6;
if(par.input), sys0.nu = 1; else sys0.nu = 0; end;
sys0.ny = 1;
sys0.type = 'Wiener';
sys0.which = '2d atan';


% Dynamics
Lorder = 6;
sys = generate_linear_system_randomly(Lorder, 0.9, 0.1);
sys = ss(sys);

sys0.ss.A = sys.a;
sys0.ss.C = sys.c;
if(par.input)
    sys0.ss.B = sys.b;
    u = randn(T,1);
else
    sys0.ss.B = [];
    u = [];
end
% 

plc = 2*rand(3,1) - 1;
% 
sys0.ss.h = @(z)   plc(1)*z.^2 + plc(2)*z + plc(3);
    

sys0.ss.Q = 0;%0.025^2*eye(sys0.nx);



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
% 
 
ytrue = sys0.ss.h(z);

snr = db2pow(10);
e = randn(T,1);
e = (e - mean(e))*std(ytrue)/std(e)/sqrt(snr);

sys0.ss.R = var(e);

y = ytrue + e;

N = T/4;
data.u = u(1:N)';
data.uv = u(N+1:2*N)';
data.ut = u(2*N+1:4*N)';

data.ytrue = ytrue(1:N)';
data.yvtrue = ytrue(N+1:2*N)';
data.yttrue = ytrue(2*N+1:4*N)';
data.z = z;
data.y = y(1:N)';
data.yv = y(N+1:2*N)';
data.yt = y(2*N+1:4*N)';
data.e = e';
data.sys = sys;
data.nonl = sys0.ss.h;

save(['databank/' 'data_' int2str(repi) '.mat'], 'data');
save(['databank/' 'sys0_' int2str(repi) '.mat'], 'sys0');
end

