clear;

c = 1;
n = 50;
ts = (1:n)';
tic
lam = 0.7;
lama =0.1;
al = -log(lama);
xi = 0.5;
bt = al/xi*sqrt(1-xi^2);
gm = -log(lam);
amg = al-gm;
phi = acos(amg/sqrt(bt^2+amg^2));


rl = exp(-al).^ts;
T = (rl*rl');

rl1 = [cos(bt*ts) sqrt((al^2+1)/bt^2)*sin(bt*ts)];
T1 = rl1*rl1';

T2 = al*sin(bt*(ts+ts'))/bt;

TE = exp(2*amg*min(ts,ts'));
T3 = cos(bt*(ts-ts')).*(TE-1)/(4*bt^2*amg);

TP = phi+bt*(ts+ts');
T4 = cos(TP)-TE.*cos(2*bt*min(ts,ts')-TP);
T4 = T4/(4*bt^2)/sqrt(bt^2+amg^2);

K = T.*(T1+T2+T3+T4);

toc
tic
hyper(1) = lama;
hyper(2) = xi;
hyper(3) = lam;
be = -log(hyper(1))/hyper(2)*sqrt(1-hyper(2)^2);
theta = acos(2*log(hyper(3)/hyper(1))/sqrt(4*be^2+4*log(hyper(3)/hyper(1))^2));
t = ts;
for dk = 1:length(t)
    for dj = 1:dk
        Pi(dk,dj) = hyper(1)^(t(dk)+t(dj))*(cos(be*t(dk))-log(hyper(1))/be*sin(be*t(dk)))*(cos(be*t(dj))-log(hyper(1))/be*sin(be*t(dj))) ...
            + hyper(1)^(t(dk)+t(dj))*sin(be*t(dk))*sin(be*t(dj))/be^2 + hyper(1)^(t(dk)+t(dj))*cos(be*(t(dk)-t(dj)))*((hyper(3)/hyper(1))^(2*min(t(dk),t(dj)))-1)/(4*be^2*log(hyper(3)/hyper(1)))...
            + hyper(1)^(t(dk)+t(dj))*(cos(theta+be*(t(dk)+t(dj))) - (hyper(3)/hyper(1))^(2*min(t(dk),t(dj)))*cos(2*be*min(t(dk),t(dj))-theta-be*(t(dk)+t(dj))) )...
            /(2*be^2*sqrt(4*be^2+4*log(hyper(3)/hyper(1))^2));
        Pi(dj,dk) = Pi(dk,dj);
    end
end
toc
K - Pi;