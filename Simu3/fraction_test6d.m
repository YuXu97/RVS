clear;

a = [-0.467, 1.12, -0.925, 0.308, -0.0364, 0.00110];
c = [-2.67, 2.96, -2.01, 0.914, -0.181, -0.0102];

sys = tf([0 a], [1 c], 1, 'Variable','z^-1');

% syms z
% 
% num = poly2sym([0 a],z);
% den = poly2sym([1 c],z);
% 
% TF = collect(num/den)

[r,p,k] = residue([0 a],[1 c]);

syms z;
[Nos,Dos] = numden(collect(r(3)/(z-p(3))+r(4)/(z-p(4))));
sys_os = tf([0,sym2poly(Nos)],sym2poly(Dos), 1, 'Variable','z^-1');


[Ndec,Ddec] = numden(collect(r(1)/(z-p(1))+r(2)/(z-p(2))+r(5)/(z-p(5))+r(6)/(z-p(6))));
sys_dec = tf([0,sym2poly(Ndec)],sym2poly(Ddec), 1, 'Variable','z^-1');


figure(1)
g = impulse(sys,0:79);
plot(0:79,g);
hold on;
g_os = impulse(sys_os,0:79,'-*');
plot(0:79,g_os);
hold on;
g_dec = impulse(sys_dec,0:79,'-o');
plot(0:79,g_dec);
grid on;
legend('g','g-os','g-dec')