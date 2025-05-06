function [return_filename] = data_generation(N, n, Morder, Lorder, snr, inputtype, noisetype, repi)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%output: Wiener type + output of a linear system
% t = (0:M-1)';
% order1 = 5;
% order2 = 5;
% order3 = 0;
%
% [system1, ~, ~] = generate_linear_system_randomly(order1, 0.95, 0.2);
% [system2, ~, ~] = generate_linear_system_randomly(order2, 0.95, 0.2);
% % numerator = [0.7568 0];
% % denominator = [1 0.5036];
% %system = tf(numerator,denominator,1);
% yL1 = lsim(system1,u,t);
% [g1,~] = impulse(system1,t);
% g1 = g1(1:n);
% yL2 = lsim(system2,u,t);
% [g2,~] = impulse(system2,t);
% g2 = g2(1:n);
% g3 = 0;
% h0 = 2;
%
% ytrue = h0 + yL1 + yL2.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output: Wiener-Hammerstein type + output of a linear system
% t = (0:M-1)';
% order1 = 5;
% order2 = 5;
% order3 = 5;
%
% [system1, ~, ~] = generate_linear_system_randomly(order1, 0.95, 0.2);
% [system2, ~, ~] = generate_linear_system_randomly(order2, 0.95, 0.2);
% [system3, ~, ~] = generate_linear_system_randomly(order3, 0.95, 0.2);
%
% yL1 = lsim(system1,u,t);
% [g1,~] = impulse(system1,t);
% g1 = g1(1:n);
%
% yL2 = lsim(system2,u,t).^2;
% yL2 = lsim(system3,yL2,t);
% [g2,~] = impulse(system2,t);
% [g3,~] = impulse(system3,t);
% g2 = g2(1:n);
% g3 = g3(1:n);
%
% h0 = 2;
%
% ytrue = h0 + yL1 + yL2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%output: Wiener system
% t = (0:M-1)';
% order = Lorder;
% a = 2*rand(Morder,1)-1;
%
% % [system, poles] = generate_linear_system_randomly(order, 0.9, 0.1);
%
% % numerator = [10 0 0];
% % denominator = [1 -1.8036 0.8338];
% % system = tf(numerator,denominator,1);
% % poles = pole(system);
%
% numerator = [2*rand-1 0];
% denominator = [1 1.9*rand-0.95]; %[1 -0.95*rand];%
% poles = -1.9*rand+0.95; %0.9*rand;
% system = tf(numerator,denominator,1);
%
% %
% % x1 = 0.9*rand;
% % x2 = 0.9*rand;
% % numerator = [2*rand-1 0 0];
% % denominator = [1 -(x1+x2) x1*x2];
% % poles = [x1 x2];
% % system = tf(numerator,denominator,1);
%
% % zeros = [0 0];
% % poles = [-inf inf];
% % poleub = 0.9; polelb = 0.1;
% % while max(abs(poles)) > poleub || min(abs(poles)) < polelb %|| max(real(poles)) > 0
% %     ac = 2*rand-1; bc = 2*rand-1; cc= 2*rand-1;
% %     poles = [-bc+sqrt(bc^2-4*ac*cc) -bc-sqrt(bc^2-4*ac*cc)]/2/ac;
% % end
% % gain = 2*rand-1;
% %
% % system = zpk(zeros,poles,gain,1);
%
%
% yL = lsim(system, u, t);
% [g,~] = impulse(system, t);
%
% h0 = 0;
%
% ytrue = h0;
%
% for i =1:Morder
%     ytrue = ytrue + a(i)*yL.^i;
% end
%
% %
% % ytrue = ytrue + a(2)*yL.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%output: Wiener-Hammerstein system
t = (0:M-1)';

a = 2*rand(Morder,1)-1;
ratioM1 = 1000;
ratioM2 = 1000;
while ratioM1 < 0.1 || ratioM1 > 10 ||...
         ratioM2 < 0.1 || ratioM2 > 10
    %     numerator1 = [2*rand-1 0];
    %     denominator1 =  [1 1.8*rand-0.9];%[1 0.8*rand+0.1];%
    %     system1 = tf(numerator1,denominator1,1);
    [system1, ~] = generate_linear_system_randomly_dd(Lorder-15, 0.9, 0.1);
    poles1 = pole(system1);
    
    
    %     numerator2 = [1 0];
    %     denominator2 = [1 -0.5]; %[1 -0.95*rand];%
    % poles2 = 0.5; %0.9*rand;
    
    % numerator2 = [2*rand-1 0];
    % denominator2 = [1 1.8*rand-0.9];% [1 0.8*rand+0.1];%
    
%     numerator2 = [0.1 0];
%     denominator2 = [1 -1.8036 0.8338];
%    
%     system2 = tf(numerator2,denominator2,1);
    [system2, ~] = generate_linear_system_randomly(Lorder, 0.9, 0.1);
    poles2 = pole(system2);
    
    yL1 = lsim(system1, u, t);
    [g1,~] = impulse(system1, t);
    
    h0 = 2*rand-1;
    z = h0;
    
    HmVar = zeros(Morder,1);
    for i =1:Morder
        z = z + a(i)*yL1.^i;
%        HmVar(i) = var(a(i)*yL1.^i);
         HmVar(i) = var(lsim(system2, a(i)*yL1.^i, t));%a(i)*yL1.^i;

    end
    
    ytrue = lsim(system2, z, t);
    [g2,~] = impulse(system2, t);
    
    ratioM1 = HmVar(2)/HmVar(1);
     ratioM2 = HmVar(3)/HmVar(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%snr = 5;%10*(0.1+0.9*rand);
switch noisetype
    case 'gauss'
        noise = randn(N,1);
    case 'uniform'
        noise = rand(N,1)-0.5;
end
noise = (noise - mean(noise))*std(ytrue)/std(noise)/sqrt(snr);
y = ytrue(end-N+1:end) + noise;
u = u(end-N+1:end);
data = [u,y];




%Wiener-Hammerstein system savings
datainfo.data = data;
datainfo.ytrue = ytrue(end-N+1:end);
datainfo.LinearSystem1 = system1;
datainfo.LinearSystemPoles1 = poles1;
datainfo.LinearSystemImpulseResponse1 = g1;
datainfo.LinearSystem2 = system2;
datainfo.LinearSystemPoles2 = poles2;
datainfo.LinearSystemImpulseResponse2 = g2;
datainfo.snr = snr;
datainfo.h0 = h0*sum(g2);
datainfo.noise = noise;
datainfo.LinearSysOrder = Lorder;
datainfo.PolyCoefficients = [h0;a];
datainfo.NoiseVariance = var(noise);
datainfo.HmVar = HmVar;


figure(1);
h = plot(0:n-1,g1(1:n)); grid on;
title(['g1-repi=' int2str(repi)]);
saveas(h,['Databank/Figs/sys1/system1_' int2str(repi) '.eps']);

figure(2);
h = plot(0:n-1,g2(1:n)); grid on;
title(['g2-repi=' int2str(repi)]);
saveas(h,['Databank/Figs/sys2/system2_' int2str(repi) '.eps']);

return_filename = ['Databank/data_N' int2str(N) '_repi=' int2str(repi) '.mat'];
save(return_filename, 'datainfo');
end


function [system, poles] = generate_linear_system_randomly(order, poleub, polelb, fs_index)
% Randomly generate a stable lienar system
% the argument 'order' is the needed order of the system
if nargin < 4 || strcmp(fs_index,'fast')
    pole_max = 2;
    pole_min = 0;
    % pole check
    value = 1;
    %     poles_real_min = -1;
    %     poles_tan_max = 10000;
    
    
    while pole_max > poleub || pole_min < polelb  || value > 1e6
        bw = inf;
        while isinf(bw)||isnan(bw) % finite bandwidth of generated system
            mc = rss(order,1,1);
            bw = bandwidth(mc);
        end
        f = 1; %bw*3*2*pi;
        md = c2d(mc,1/f,'zoh');
        md.d = 2*rand-1;
        system = idpoly(md);
        
        poles = pole(system);
        pole_max = max(abs(poles));
        pole_min = min(abs(poles));
        
        value = max(abs([system.f system.b]));
    end
    
    
    
    
    % five poles with the largest modulus is real and around 0.75
    % pole_max = 100;
    % tmp = 100;
    % while tmp > 1
    %     while pole_max > poleub || pole_min < polelb  || value > 1e6 || pole_max > 0.4
    %         bw = inf;
    %         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
    %             mc = rss(order-5,1,1);
    %             bw = bandwidth(mc);
    %         end
    %         f = 1;%bw*3*2*pi;
    %         md = c2d(mc,1/f,'zoh');
    %         md.d = 0;
    %         system_t = idpoly(md);
    %         poles = pole(system_t);
    %         pole_max = max(abs(poles));
    %         pole_min = min(abs(poles));
    %         value = max(abs([system_t.f system_t.b]));
    %     end
    %     zeros_zpk = [zero(system_t);rand;rand(5,1)];%make sure there is zeros delay
    %     poles_zpk = [poles; 0.2*rand(5,1)+0.7];
    %     system = zpk(zeros_zpk,poles_zpk,2*rand-1,1);
    %     [g,~] = impulse(system, 0:10);
    %
    %     tmp = max(abs(g));
    % end
    %
    %
    
    
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

function [system, poles] = generate_linear_system_randomly_dd(order, poleub, polelb, fs_index)
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
    %
    %    while pole_max > poleub || pole_min < polelb  || value > 1e6
    %        bw = inf;
    %         while isinf(bw)||isnan(bw) % finite bandwidth of generated system
    %             mc = rss(order,1,1);
    %             bw = bandwidth(mc);
    %         end
    %         f = 1; %bw*3*2*pi;
    %         md = c2d(mc,1/f,'zoh');
    %         md.d = 2*rand-1;
    %         system = idpoly(md);
    %
    %         poles = pole(system);
    %         pole_max = max(abs(poles));
    %         pole_min = min(abs(poles));
    %
    %         value = max(abs([system.f system.b]));
    %    end
    
    
    
    
    % five poles with the largest modulus is real and around 0.75
    pole_max = 100;
    
    
    while pole_max > poleub || pole_min < polelb  || value > 1e6 || pole_max > 0.5
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
    zeros_zpk = [zero(system_t);2*rand(5,1)-1];%make sure there is a delay
    poles_zpk = [poles; 0.1*rand(5,1)+0.7];
    system = zpk(zeros_zpk,poles_zpk,2*rand-1,1);
    %[g,~] = impulse(system, 0:10);
    
    
    
    %
    %
    
    
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

