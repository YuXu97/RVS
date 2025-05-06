function [opk_gen] = create_opk(kernel_gen,inp_gen,Ts)
%% creat generatrors for output kernel

% Inputs
% kernel_gen: the generators of the kernel
% inp_gen: the generators of the input
% T: sampling time, set to 1 if not specified.
% Nest: the number of time instants estimated for the impulse response
%
% Outputs
% opk_gen: the generators of the output kernel


assert(mod(size(inp_gen,2),2) == 0,'the number of generators of the input should be even');
assert(mod(size(kernel_gen,2),2) == 0,'the number of generators of the kernel should be even');
assert(size(inp_gen,1) == size(kernel_gen,1),'the number of points in the generator of the input should be same as the number of data');
assert(size(kernel_gen,1) == size(kernel_gen,1),'the number of points in the generator of the kernel should be same as the number of data');

N = size(kernel_gen,1);
r = size(inp_gen,2)/2;
p = size(kernel_gen,2)/2;


% retrive the generator representation of the kernel and the input
mu = kernel_gen(:,1:p);
nu = kernel_gen(:,p+1:end);

pai = inp_gen(:,1:r);
rho = inp_gen(:,r+1:end);

if nargin <= 2
    Ts=1;
end

bU = zeros(N,p+r);
bV = zeros(N,p+r);


murho = zeros(N,p*r);
for j=1:r
    for i=1:p
        murho(:,i+(j-1)*p) = mu(:,i).*rho(:,j);%murho(:,i+(j-1)*p):mu_i*rho_j
    end
end
% f_{ij}^1(t) = \int_0^t \mu_i(a)\rho_j(a)da
csum_murho = cumsum(murho,1);


nurho = zeros(N,p*r);
for j=1:r
    for i=1:p
        nurho(:,i+(j-1)*p) = nu(:,i).*rho(:,j);%nurho(:,i+(j-1)*p):nu_i*rho_j
    end
end
% f_{ij}^2(s) = \int_0^s \nu_i(a)\rho_j(a)da
csum_nurho = cumsum(nurho,1);


%\bar\mu_i(t), i=1,...,p
% \sum_{j=1}^r \pi_j(t)f_{i,j}^1(t)
for i=1:p
    tmpj=0;
    for j=1:r
        tmpj = tmpj + pai(:,j).*csum_murho(1:N,i+(j-1)*p)*Ts;
        % csump2(1:N,i+(j-1)*p) : f_ij(t), f_ij(s)
    end
    
    bU(:,i)=tmpj;  %the U matrix with 1st to the pth columns
    
end

% \bar\nu_i(s), i=1,...,p
% \bar\nu_i(s) = \sum_{j=1}^r \pi_j(s) f_{i,j}^2(s)
for i=1:p
    tmpj=0;
    for j=1:r
        tmpj = tmpj + pai(:,j).*csum_nurho(1:N,i+(j-1)*p)*Ts;
        % csump1(1:N,i+(j-1)*p)
        
    end
    bV(:,i)=tmpj;    %the V matrix with 1st to the pth columns
    
end


% \bar\mu_i(t), i=p+1,...,p+r,
% \pi_i(t)
for j=1:r
    bU(:,j+p) = pai(:,j); %the U matrix with (p+1)th to the (p+r)th columns
end

% \bar\nu_i(s), i=p+1,...,p+r,
% \bar\ell_j(s) = \sum_{i=1}^p \bar\nu_i(s)(-f_{ij}^1(s))
for j=1:r
    
    tmpi=0;
    for i=1:p
        tmpi = tmpi - bV(:,i).*csum_murho(1:N,i+(j-1)*p)*Ts;
        
    end
    
    bV(:,j+p) = tmpi; %the V matrix with (p+1)th to the (p+r)th columns
    %bb(:,j) = tmpi;
end


%\bar\nu_i(s), i=p+1,...,p+r
%\bar\rho_i(s), i=1,...,r
for i=1:r
    
    tmpl = 0;
    for l=1:r
        
        tmpj = 0;
        for j=1:p
            tmpj = tmpj + ...
                cumsum([0; csum_nurho(1:N-1,j+(l-1)*p)].*murho(1:N,j+(i-1)*p) +...
                nurho(1:N,j+(l-1)*p).*murho(1:N,j+(i-1)*p))+...
                cumsum([0; 0; csum_nurho(1:N-2,j+(i-1)*p)].*murho(1:N,j+(l-1)*p) +...
                murho(1:N,j+(l-1)*p).*[0; nurho(1:N-1,j+(i-1)*p)]);
        end 
        tmpl = tmpl + pai(:,l).*tmpj*Ts^2;
    end
    bV(:,i+p) = bV(:,i+p) + tmpl;
    
end

opk_gen = [bU, bV];

%Psi = tril(bU*bV')+triu(bV*bU',1); for test
    
end