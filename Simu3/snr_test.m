
a = [-0.467, 1.12, -0.925, 0.308, -0.0364, 0.00110];
c = [-2.67, 2.96, -2.01, 0.914, -0.181, -0.0102];

sys = tf([0 a], [1 c], 1, 'Variable','z^-1');
sys = ss(sys);



snr = zeros(40,1);
maxr = 50000;
for repi = 1:maxr
    u = randn(10000,1);
    z = lsim(sys, u, 1:10000); 
    h = @(z)(2*max(min(z,0.5),-0.5));
%     e = 0.1*randn(1000,1);
    y = h(z);
    snr(repi) = pow2db(var(y)/0.01); 
end

plot(1:maxr,snr)
mean(snr)