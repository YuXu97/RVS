clear;

xx = (-4:1e-2 :4)';

figure(1)

h = @(z)nonmon(z,0.5,0.3);

plot(xx,h(xx))
hold on;


order = 13;
pl = polyfit(xx,h(xx),order);

poly = @(z) polyval(pl,z);

plot(xx, poly(xx));
grid on;

T = 500;

u = 0.5*randn(T,1);



c= [0.368, 0.888, 0.524, 0.555];
a = [1, 0.1, -0.49, 0.01];

sys = tf([0 c], [1 a], 1, 'Variable','z^-1');


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

M = 7;
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