clear;

xx = (-2:1e-2 :2)';

figure(1)
h = @(z) (2*max(min(z,0.5),-0.5));


plot(xx,h(xx))
hold on;


order = 11;
pl = polyfit(xx,h(xx),order);

poly = @(z) polyval(pl,z);

plot(xx, poly(xx));
grid on;

T = 500;

u = randn(T,1);

a = [-0.467, 1.12, -0.925, 0.308, -0.0364, 0.00110];
c = [-2.67, 2.96, -2.01, 0.914, -0.181, -0.0102];

sys = tf([0 a], [1 c], 1, 'Variable','z^-1');


g = impulse(sys,0:T-1);
n = 150;
g = g(1:n);
z = lsim(sys, u, 1:T); 
zt = CalculatePsi(u,n)*g;
d = z(n:end)-zt;

figure(2)
plot(z)
(h(z) - poly(z))./h(z);


id = 1:2:order+1;

coef = pl(flip(id));

M = 5;
pll = zeros(T,1);
for i = 1:M
    pll = pll + coef(i)*z.^(2*i-1);
end

pll2 = zeros(T,1);
TMP = z;
DP = z.^2;
for i = 1:M
    pll2 = pll2 + coef(i)*TMP;
    TMP = TMP.*DP;
end

pll - pll2;


gof(h(z),pll)
save('g.mat','g')
save('coef.mat','coef')