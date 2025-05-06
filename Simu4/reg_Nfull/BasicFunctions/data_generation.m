function [return_filename] = data_generation(N, n, Morder, Lorder, snr, inputtype, poleub, noisetype, noisevar, repi)
warning('off');
M = 1000 + N;

switch inputtype{1}

    case 'white'
        fd_range = inputtype{2};
        fb_low = fd_range(1);
        fd_upper = fd_range(2);
        u = idinput(M,'rgs',[fb_low,fd_upper],[]);

    case 'filtered'
        Syst_filter = generate_linear_system_randomly(30,0.95,0.1);
        u = randn(M,1);
        u = sim(Syst_filter,u);
        
    case 'uniform'
        u = rand(M,1);%2*rand(M,1)-1;
end


%%output: Wiener system
t = (0:M-1)';
order = Lorder;
% 
a = ones(Morder,1);
%      a = 5*(2*rand(Morder,1)-1);
%  a = rand(Morder,1);
% ratioM = 1000;
% while ratioM < 0.5 || ratioM > 2
% [system, ~] = generate_linear_system_randomly(order, 0.95, 0.1);
% % 
% numerator = [2*rand-1 0];
% denominator = [1 -0.95*rand];%[1 -0.95*rand];
% numerator = [0.3 0];
% denominator = [1 -1.1232 0.8099];%[1 -0.95*rand];
numerator = [0.06 0];
denominator = [1 -1.8036 0.8338];
% 
% numerator = [0.1 0];
% denominator = [1 -1.5 0.81];


% numerator = [0.1 0];
% denominator = [1 -0.3 0.52];

% a1 = 0.95*(2*rand-1); b1 = 0.95*(2*rand-1);
% numerator = [2*rand-1 0];
% denominator = [1 -(a1+b1) a1*b1];
% ratio = -1;
% while ratio < 0
%     rho = rand; theta = pi/2*rand;
%     a1 = rho*cos(theta)+1i*rho*sin(theta);
%     b1 = rho*cos(theta)-1i*rho*sin(theta);
%     ratio = cos(theta);
% end
% numerator = [2*rand-1 0];
% denominator = [1 -(a1+b1) a1*b1];


system = tf(numerator,denominator,1);
poles = pole(system);



% numerator1 = [10 0];
% denominator1 = [1 -1.8036 0.8338];
% 
% system1 = tf(numerator1,denominator1,1);
% 
% numerator2 = [2*rand-1 0];
% denominator2 = [1 0.95*rand]; %[1 -0.95*rand];%t
% % poles = -1.9*rand+0.95; %0.9*rand;
% system2 = tf(numerator2,denominator2,1);






yL = lsim(system, u, t); 
[g,~] = impulse(system, t);

h0 = 0;

ytrue = h0;
VarM = zeros(Morder,1);
for i =1:Morder
    ytrue = ytrue + a(i)*yL.^i;
    VarM(i) = var(a(i)*yL.^i);
end
% ratioM = VarM(2)/VarM(1);
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%snr = 5;%10*(0.1+0.9*rand);
switch noisetype
    case 'gauss'
        noise = randn(N,1);
    case 'uniform'
        noise = rand(N,1)-0.5;
end
noise = (noise - mean(noise))*std(ytrue)/std(noise)/sqrt(snr);
% noise = noise/std(noise)*sqrt(noisevar);


y = ytrue(end-N+1:end) + noise;
u = u(end-N+1:end);
data = [u,y];

% datainfo.data = data;
% datainfo.ytrue = ytrue(end-N+1:end);
% datainfo.LinearSystemImpulseResponse1 = g1;
% datainfo.LinearSystemImpulseResponse2 = g2;
% datainfo.LinearSystemImpulseResponse3 = g3;
% datainfo.snr = snr;
% datainfo.noise = noise;
% datainfo.LinearSysOrder1 = order1;
% datainfo.LinearSysOrder2 = order2;
% datainfo.LinearSysOrder3 = order3;
% datainfo.NoiseVariance = var(noise);

datainfo.data = data;
datainfo.ytrue = ytrue(end-N+1:end);
datainfo.LinearSystem = system;
datainfo.LinearSystemPoles = poles;
datainfo.LinearSystemImpulseResponse = g;
datainfo.LinearSystemOutput = yL(end-N+1:end);
datainfo.h0 = h0;
datainfo.snr = snr;
datainfo.noise = noise;
datainfo.LinearSysOrder = order;
datainfo.PolyCoefficients = [h0;a];
datainfo.NoiseVariance = var(noise);
datainfo.VarianceOfEachOrder = VarM;

figure(1);
h = plot(0:n-1,g(1:n)); grid on;
title(['repi=' int2str(repi)]);
saveas(h,['Databank/Figs/system' int2str(repi) '.eps']);

return_filename = ['Databank/data_N' int2str(N) '_repi=' int2str(repi) '.mat'];
save(return_filename, 'datainfo');
end


function [system, poles] = generate_linear_system_randomly(order, poleub, polelb, fs_index)
% Randomly generate a stable lienar system
% the argument 'order' is the needed order of the system
if nargin < 4 || strcmp(fs_index,'fast')
%     pole_max = 2;
%     pole_min = 0;
%     % pole check
%     value = 1;
% %     poles_real_min = -1;
% %     poles_tan_max = 10000;
% 
% % ratio = 0;
%    while pole_max > poleub || pole_min < polelb  || value > 1e6 %|| ratio < 0.8 
%        bw = inf;
%         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
%             mc = rss(order,1,1);
%             bw = bandwidth(mc);
%         end
%         f = 1; %bw*3*2*pi;
%         md = c2d(mc,1/f,'zoh');
%         md.d = 0;
%         system = idpoly(md);
%         
%         poles = pole(system);
%         pole_max = max(abs(poles));
%         pole_min = min(abs(poles));
%                 
%         value = max(abs([system.f system.b]));
% %         ratio = cos(angle(poles(1)));
% 
%    end
% 

% %poles that are close to the real line
% pole_max = 2;
% pole_min = 0;
% poles_imag_max = 2;
% 
% while pole_max > poleub || pole_min < polelb  || value > 1e6 || poles_imag_max > 0.1
%     bw = inf;
%     while isinf(bw)||isnan(bw) % finite bandwidth of generated system
%         mc = rss(order,1,1);
%         bw = bandwidth(mc);
%     end
%     f = 1; %bw*3*2*pi;
%     md = c2d(mc,1/f,'zoh');
%     md.d = 0;
%     system = idpoly(md);
%     
%     poles = pole(system);
%     pole_max = max(abs(poles));
%     pole_min = min(abs(poles));
%     poles_imag_max = max(imag(poles));
%     
%     value = max(abs([system.f system.b]));
% end
% system = zpk([zero(system);0],pole(system),2*rand-1,1);


   
% five poles with the largest modulus is real and around 0.75
pole_max = 100;

while pole_max > poleub || pole_min < polelb  || value > 1e6 || pole_max > 0.45
    bw = inf;
    while isinf(bw)||isnan(bw) % finite bandwidth of generated system
        mc = rss(order-5,1,1);
        bw = bandwidth(mc);
    end
    f = 1;%bw*3*2*pi;
    md = c2d(mc,1/f,'zoh');
    md.d = 0;
    system_t = idpoly(md);
    poles = pole(system_t);
    pole_max = max(abs(poles));
    pole_min = min(abs(poles));
    value = max(abs([system_t.f system_t.b]));
end
zeros_zpk = [zero(system_t);rand;rand(5,1)];%make sure there is zeros delay
poles_zpk = [poles; 0.25*rand(5,1)+0.7];
system = zpk(zeros_zpk,poles_zpk,2*rand-1,1);


   
   
   
%    num_con = 0;
% 
%     while pole_max > poleub || pole_min < polelb  || value > 1e6 || num_con < ceil(order*4/5)
%         bw = inf;
%         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
%             mc = rss(order,1,1);
%             bw = bandwidth(mc);
%         end
%         f = 1;%bw*3*2*pi;
%         md = c2d(mc,1/f,'zoh');
%         md.d = 2*rand-1;
%         system = idpoly(md);
%         
%         poles = pole(system);
%         pole_max = max(abs(poles));
%         pole_min = min(abs(poles));
%         
%         con1 = find(real(poles)<0);
%         con2 = con1;
%         num_con = length(intersect(con1,con2));
%         
%         value = max(abs([system.f system.b]));
%     end


%    num_con = 0;
% 
%     while pole_max > poleub || pole_min < polelb  || value > 1e6 || num_con < ceil(order)
%         bw = inf;
%         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
%             mc = rss(order,1,1);
%             bw = bandwidth(mc);
%         end
%         f = 1;%bw*3*2*pi;
%         md = c2d(mc,1/f,'zoh');
%         md.d = rand;
%         system = idpoly(md);
%         
%         poles = pole(system);
%         pole_max = max(abs(poles));
%         pole_min = min(abs(poles));
%         
%         pole_real = real(poles);
%         pole_imag = imag(poles);
%         con1 = find(pole_real>0);
%         con2 = find(pole_imag==0)
%         num_con = length(intersect(con1,con2));
%         
%         value = max(abs([system.f system.b]));
%     end
    
    
    poles = pole(system);
    num_poles = length(pole(system));
    num_zeros = length(zero(system));
    pz = [num_poles num_zeros];


elseif strcmp(fs_index,'slow')

    pole_max = 0.5; % pole check
    value = 1;
    while pole_max < poleub || pole_max > 1 - eps || value > 1e6
        bw = inf;
        while isinf(bw)||isnan(bw) % finite bandwidth of generated system
            mc = rss(order,1,1);
            bw = bandwidth(mc);
        end
        f = bw*3*2*pi;
        md = c2d(mc,1/f,'zoh');
        md.d = 0;
        system = idpoly(md);
        pole_max = max(abs(pole(system)));
        value = max(abs([system.f system.b]));
    end
end

end
