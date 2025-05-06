clear;

c= [0.368, 0.888, 0.524, 0.555];
a = [1, 0.1, -0.49, 0.01];

sys = tf([0 c], [1 a], 1, 'Variable','z^-1');


[r,p,k] = residue([0 c],[1 a]);

syms z;
[Nos,Dos] = numden(collect(r(1)/(z-p(1))+r(2)/(z-p(2))));
sys_os = tf([0,sym2poly(Nos)],sym2poly(Dos), 1, 'Variable','z^-1');
 

[Ndec,Ddec] = numden(collect(r(3)/(z-p(3))+r(4)/(z-p(4))));
sys_dec = tf([0,sym2poly(Ndec)],sym2poly(Ddec), 1, 'Variable','z^-1');
% 
% 
figure(1)
g = impulse(sys,0:150);
plot(0:150,g);
hold on;
g_os = impulse(sys_os,0:150);
plot(0:150,g_os,'-*');
hold on;
g_dec = impulse(sys_dec,0:150);
plot(0:150,g_dec,'-o');
grid on;
legend('g','g-os','g-dec')