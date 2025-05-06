function [system, poles] = generate_linear_system_randomly(order, poleub, polelb, fs_index)
% Randomly generate a stable lienar system
% the argument 'order' is the needed order of the system
if nargin < 4 || strcmp(fs_index,'fast')
    pole_max = 2;
    pole_min = 0;
    % pole check
    value = 1;

   while pole_max > poleub || pole_min < polelb  || value > 1e6 %|| ratio < 0.8 
       bw = inf;
        while isinf(bw)||isnan(bw) % finite bandwidth of generated system
            mc = rss(order,1,1);
            bw = bandwidth(mc);
        end
        f = 1; %bw*3*2*pi;
        md = c2d(mc,1/f,'zoh');
        md.d = 0;
        system = idpoly(md);
        
        poles = pole(system);
        pole_max = max(abs(poles));
        pole_min = min(abs(poles));
                
        value = max(abs([system.f system.b]));
   end


    
    
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
