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
        
    case 'filtered_multisine'
        Period = M;
        NumPeriod = 1;
        Range = [-1 1];
        band = [0 1];
        NumSin = 100;
        NumTrials = 20;
        GridSkip = 2;
        SinData = [NumSin,NumTrials,GridSkip];
        u = idinput([Period 1 NumPeriod],'sine',band,Range,SinData);
        Syst_filter = generate_linear_system_randomly(2,0.95,0.1);
        u = sim(Syst_filter,u);
        u = u/rms(u);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output: Wiener-Hammerstein type + output of a linear system
t = (0:M-1)';

numerator1 = [0.7568 0];
denominator1 = [1 -1.812 0.8578];

system1 = tf(numerator1,denominator1,1);
poles1 = pole(system1);

numerator2 = [1.063 0];
denominator2 = [1 -1.706 0.7491];

system2 = tf(numerator2,denominator2,1);
poles2 = pole(system2);

% system3 = generate_linear_system_randomly(15,0.9,0.1);
system3 = tf(1.5*numerator1,denominator1,1);
poles3 = pole(system3);

yL1 = lsim(system1,u,t);
[g1,~] = impulse(system1,t);
g1 = g1(1:n);

yL2 = lsim(system2,u,t).^2;
yL2 = lsim(system3,yL2,t);
[g2,~] = impulse(system2,t);
[g3,~] = impulse(system3,t);
g2 = g2(1:n);
g3 = g3(1:n);

h0 = 2;

ytrue = h0 + yL1 + yL2;
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


%Wiener-Hammerstein type + linear system savings

H = zeros(n,n);

for t1 = 1:n
    for t2 = 1:n
        for  sig = 1:n
            if t1 >= sig && t2 >= sig
                H(t1,t2) = H(t1,t2) + g3(sig)*g2(t1-sig+1)*g2(t2-sig+1);
            end
        end
    end
end
save('H.mat','H');

datainfo.SecondOrderVolterraKernel = H;
datainfo.data = data;
datainfo.ytrue = ytrue(end-N+1:end);
datainfo.LinearSystem1 = system1;
datainfo.LinearSystemPoles1 = poles1;
datainfo.LinearSystemImpulseResponse1 = g1;
datainfo.LinearSystem2 = system2;
datainfo.LinearSystemPoles2 = poles2;
datainfo.LinearSystemImpulseResponse2 = g2;
datainfo.LinearSystem3 = system3;
datainfo.LinearSystemPoles3 = poles3;
datainfo.LinearSystemImpulseResponse3 = g3;
datainfo.snr = snr;
datainfo.h0 = h0;
datainfo.noise = noise;
datainfo.LinearSysOrder = Lorder;
datainfo.NoiseVariance = var(noise);



% h = surf(H);
% 
% figure(1);
% h = plot(0:n-1,g1(1:n)); grid on;
% title(['g1-repi=' int2str(repi)]);
% saveas(h,['Databank/Figs/sys1/system1_' int2str(repi) '.eps']);
% 
% figure(2);
% h = plot(0:n-1,g2(1:n)); grid on;
% title(['g2-repi=' int2str(repi)]);
% saveas(h,['Databank/Figs/sys2/system2_' int2str(repi) '.eps']);

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
    poles_real_min = -1;
    poles_tan_max = 10000;


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


    
   
% % five poles with the largest modulus is real and around 0.75
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
% %    
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
