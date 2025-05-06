
addpath([pwd '/BasicFunctions']);
addpath([pwd '/conv2fft']);
addpath([pwd '/TuningMethods']);

N = 379;
data = datainfo.data(1:N,:);


n = 80;
u = data(:,1);
Psi = CalculatePsi(u, n);
g1 = datainfo.LinearSystemImpulseResponse1;
g2 = datainfo.LinearSystemImpulseResponse2;
g3 = datainfo.LinearSystemImpulseResponse3;
% Psi = Psi(n:end,:);
yconv = conv(g3,(Psi*g2).^2,'full');
yconv = yconv(n+1:end-n+1);
yconvt = conv((Psi*g2).^2, g3,'valid');
yconvt = yconvt(2:end);


numerator1 = [0.7568 0];
denominator1 = [1 -1.812 0.8578];
system1 = tf(numerator1,denominator1,1);

numerator2 = [1.063 0];
denominator2 = [1 -1.706 0.7491];
system2 = tf(numerator2,denominator2,1);

system3 = tf(1.5*numerator1,denominator1,1);

t = 0:N-1;
yL1 = lsim(system1,u,t);
yL2 = lsim(system2,u,t).^2;
yL2 = lsim(system3,yL2,t);


figure(1)
plot(yconvt); hold on;plot(yconv); hold on;; plot(yL2(2*n:N)); legend('yconvt','yconv','yL2')

% figure(2)
% plot(Psi(n+1:end,:)*g1); hold on; plot(yL1(2*n:N)); legend('yconv','yL2')

