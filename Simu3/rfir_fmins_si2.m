function  [M, estinfo] = rfir_fmins_si2(varargin)
%new_rfir	Construct an FIR model using regularization method.
%
%   M = new_rfir(z,nb,kernel) or
%   M = new_rfir(z,[nb,nk],kernel) or
%   M = new_rfir(z,{nb,nk},kernel)
%
%   estimate an FIR model represented by:
%   y(t) = B(q) u(t-nk) +  e(t)
%   where:
%       nb = order of B polynomial + 1
%       nk = input delay (in number of samples)
%       Ny = number of output
%       Nu = number of inputs
%
%   The estimated model is delivered as an IDPOLY object. Type "help
%   idpoly" for more information on IDPOLY objects.
%
%   Output:
%       M : IDPOLY model containing estimated values for B
%       polynomials along with their structure information. When Ny = 1,
%       M is a multi-input idpoly model. When Ny > 1, Jl is temporarily
%       a Ny by 1 cell array with each element as a multi-input idpoly model.
%       In the future version, Jl will be a multi-input-multi-ouput idpoly model.
%
%   Inputs:
%       z : The estimation data as an IDDATA. Use IDDATA object for input-output
%       time domain signals. Type "help iddata" for more information.
%
%       nb and nk: Orders and delays of the FIR model. nb must be specified
%       but nk can be omitted and the default value of nk is 1. Both nb and
%       nk can be scalar, which assumes all FIR models will have the same
%       order and input delay. nb can also be a Ny by 1 column vector, each
%       element of which indicates the order of all FIR models associated
%       with that ouput. nk can also be a Ny by Nu matrix, through which the
%       delays of each FIR model can be set individually.
%
%       kernel: the structure of the regularization matrix. The following
%       kernels are supported:
%            CS: qubic spline kernel
%            SE: squared exponential kernel
%            SS: stable spline kernel
%            HF: high frequency stable spline kernel
%            DI: diagonal kernel
%            TC: tuned and correlated kernel
%            DC: diagonal and correlated kernel
%       kernel can be a char (e.g., kernel = 'DC'), which assumes all FIR models
%       use the same kernel. kernel can also be a Ny by Nu cell array,
%       through which the kernels of each FIR model can be set individually.
%       More information about the structure of the regularization matrix and
%       the regulariation method can be found in literature.
%
%
%  Author: Tianshi Chen
% check No. of inputs
narginchk(3,3); % data, [nb nk], kernel

% front matter
% decomopse the inputs and check applicability
for j = 1:length(varargin)
    if isa(varargin{j},'iddata')
        data = varargin{j};
    else if isscalar(varargin{j}) || isnumeric(varargin{j})
            nbnk = varargin{j};
        else if ischar(varargin{j})
                kernel = upper(varargin{j});
            else if iscell(varargin{j})
                    tmp = varargin{j};
                    if isnumeric(tmp{1})
                        nbnk = tmp;
                    else
                        kernel = upper(tmp);
                    end
                end
            end
        end
    end
end

if realdata(data)
    [Nd, Ny, Nu, Ne] = size(data);
    dom = pvget(data,'Domain');
else
    error('Error in the data: only real-valued signals supported.');
end



if isscalar(nbnk)
    nb = nbnk*ones(Ny,Nu); nk = ones(Ny,Nu);
end

if isnumeric(nbnk) && ~isscalar(nbnk)
    if isequal(size(nbnk),[1 2])
        nb = nbnk(1)*ones(Ny,Nu); nk = nbnk(2)*ones(Ny,Nu);
    else if isequal(size(nbnk),[Ny Nu])
            nb = nbnk; nk = ones(Ny,Nu);
        else
            error('Error in the order (nb) and delay (nk) parameter matrix: Dimension mismatch with that of the MIMO system.')
        end
    end
end

if iscell(nbnk)
    if isequal(size(nbnk{1}),[Ny Nu]) && isequal(size(nbnk{2}),[Ny Nu])
        nb = nbnk{1}; nk = nbnk{2};
    else
        error('Error in the order (nb) and delay (nk) parameter matrix: Dimension mismatch with that of the MIMO system.')
    end
end


if ischar(kernel)
    tmp = kernel;
    kernel = cell(Ny,Nu);
    for i = 1:Nu
        for j = 1:Ny
            kernel{j,i} = tmp;
        end
    end
else if iscell(kernel)
        if ~isequal(size(kernel),[Ny Nu])
            error('Error in the kernel allocation matrix: Dimension mismatch with that of the MIMO system.')
        end
    end
end

indth = zeros(Ny,Nu+1);
for j = 1:Ny
    ind = cumsum(nb(j,:)+1-nk(j,:));  % ind: index of the estimated FIR parameters for each of Nu blocks
    indth(j,:) = [0 ind];
    
    for i = 1:Nu
        if ~strcmp(kernel{j,i},{'CS' 'SE' 'DI' 'SS' 'HF' 'TC' 'DC' 'NEW1' 'NEW11' 'NEW2' 'SI2' 'ML1' 'ML2' 'ML3' 'SSP'})
            error(['Error in output channel No.' num2str(j) ': No such kernel ' kernel{j,i} ' supported.']);
        end
    end
end



% main body

th = cell(Ny,1);
hp = cell(Ny,1);
hpini = cell(Ny,1);
obj = zeros(Ny,1);
exflg = zeros(Ny,1);
sigma = zeros(Ny,1);
indhp = zeros(Ny,Nu+1);
% gradient = cell(Ny,1);
% hessian = cell(Ny,1);
% otpt = cell(Ny,1);
condNb = zeros(Ny,1);

% alg = {'trust-region-reflective','active-set', 'interior-point','sqp'};
% options = optimset('GradObj','on','Hessian','on','Hessian','user-supplied','Algorithm',alg{1},'Display','off');

for j = 1:Ny
    z = data(:,j,:);
    [hpini{j}, Rt, Ney, sigma(j), indhp(j,:),lb,ub] = ini_rfir_miso(z,nb(j,:),nk(j,:),kernel(j,:));
    ff = @(x)nglglklhd_simp(x,Rt, Ney, kernel(j,:),sigma(j),indth(j,:), indhp(j,:),lb,ub);
    [hp{j}, obj(j), exflg(j)] = fminsearch(ff,hpini{j});
    %      ff = @(x)nglglklhd_simp(x,Rt, Ney, kernel(j,:),sigma(j),indth(j,:), indhp(j,:));
    %      [hp{j} obj(j) exflg(j) otpt{j}] = fminsearch(ff,hpini{j});
    [~, condNb(j), th{j}] = nglglklhd_simp(hp{j},Rt, Ney, kernel(j,:),sigma(j),indth(j,:), indhp(j,:),lb,ub);
end


% output the estimated model and associated structure information
if Ny == 1
    M =  idpoly();
else
    M = idarx();
end

for j = 1:Ny
    
    thtmp = th{j};
    mtmp = idpoly();
    
    for i = 1:Nu
        md = idpoly(1,[zeros(1,nk(j,i)) thtmp(indth(j,i)+1:indth(j,i+1))']);
        mtmp = [mtmp md];
    end
    
    mtmp.NoiseVariance =  sigma(j);
    estinfo = mtmp.estimationinfo;
    estinfo.Status =  ['estimated ' 'model'];
    estinfo.Method =  'new_rfir';
    estinfo.DataLength =  Nd;
    estinfo.DataDomain =  'Time';
    estinfo.InitialState ='Zero';
    estinfo.IniHyperparameter = hpini{j}';
    estinfo.Hyperparameter = hp{j}';
    estinfo.Hyperpartition = indhp(j,:)';
    estinfo.Kernel = [kernel{j,:}];
    %     estinfo.Optimizationinfo = struct('Minimizer','fmincon', 'options', options, 'Gradient',gradient{j},'Hessian', hessian{j},'Exitflag', exflg(j),'Optinfo', otpt{j});
    estinfo.obj = obj(j);
    estinfo.cond =  condNb;
    try
        mtmp.est = estinfo;
    catch
    end
    if Ny == 1
        M = mtmp;
    else
        if j == 1
            M = idarx(mtmp);
        else
            M = [M;idarx(mtmp)];
        end
    end
end



function [hpini, Rt, Ney, sigma, ind2, lb, ub] = ini_rfir_miso(z,nb,nk,kernel,options)

Nu = size(z,3);

% construct the regression matrix Phi for sys-id
[Rt, Ney, sigma] = qrfactor_Phi(z,nb,nk);


% Initialize the prior hyperparameter
ind2 = zeros(1,Nu);
hpini = [];
lb = [];
ub = [];
tol = sqrt(eps);
for ni = 1:Nu
    
    if strcmp(kernel{ni}, 'CS')
        lbi = -inf;
        ubi = +inf;
    end
    if strcmp(kernel{ni}, 'SE')
        lbi = [tol  -inf]';
        ubi = [+inf +inf]';
    end
    if sum(strcmp(kernel{ni}, {'SS' 'HF' 'DI' 'TC' 'NEW1' 'NEW11'}))
        lbi = [tol -inf]';
        ubi = [1-tol +inf]';
    end
    if strcmp(kernel{ni}, 'SSP')
        lbi = [ tol tol 0 -inf]';% decay, modulus of the pole of the parametric part, angle, scale
        ubi = [ 1-tol 1-tol 2*pi +inf]';
    end
    if strcmp(kernel{ni}, 'DC')
        lbi = [tol -1+tol -inf]';
        ubi = [1-tol 1-tol +inf]';
    end
    if strcmp(kernel{ni}, 'NEW2')
        lbi = [tol tol -inf]';
        ubi = [1-tol 1-tol +inf]';
    end
    if strcmp(kernel{ni}, 'SI2')
        lbi = [ tol tol tol -inf]';
        ubi = [ 1-tol 1-tol 1-tol +inf]';
    end
    if strcmp(kernel{ni}, 'ML1')
        lbi = [ tol -inf  -inf]'; %lambda alpha(freq) c
        ubi = [ 1-tol +inf +inf]';
    end
    
    if strcmp(kernel{ni}, 'ML2')
        lbi = [ tol -inf -1+tol -inf]'; %lambda omega(freq,decay) rho c
        ubi = [ 1-tol +inf 1-tol +inf]';
    end
    
    if strcmp(kernel{ni}, 'ML3')
        lbi = [ tol -inf -inf -inf]'; % lambda omega(freq,decay) alpha(freq,corr)  c
        ubi = [ 1-tol +inf +inf +inf]';
    end
    
    
    if Nu == 1
        Rti = Rt;
    else
        Rti = qrfactor_Phi(z(:,:,ni),nb(ni),nk(ni));
    end
    
    if ~sum(strcmp(kernel{ni},{'NEW2' 'SI2' 'ML1' 'ML2' 'ML3' 'SSP'}))
        %     beta =[0.5:0.1:0.9 0.99]';
        %     hpr = [-6 -4 -2 0:1:5];
        beta =[0.5:0.01:0.9 0.95]';
        hpr = -6:0.5:15;
        obj = zeros(size(beta,1),size(hpr,2));
        for nj = 1:size(beta,1)
            for nm = 1:size(hpr,2)
                hpstart = [beta(nj)*ones(1,~strcmp(kernel{ni}, 'CS')) 0.99*ones(1,strcmp(kernel{ni}, 'DC')) hpr(nm)]';
                obj(nj,nm)
        end = nglglklhd_simp(hpstart,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)],lbi,ubi);
            end
        
        [~, ind_bbeta] = min(obj(:));
        [indr, indc] = ind2sub([nj nm],ind_bbeta);
        hptmp = [beta(indr)*ones(1,~strcmp(kernel{ni}, 'CS')) 0.99*ones(1,strcmp(kernel{ni}, 'DC')) hpr(indc)]';
    end
    
    %     if strcmp(kernel{ni},'DC')
    %         beta =[0.5:0.01:0.9 0.95];
    %         rho =[-0.99 -0.9:0.1:0.9 0.99];
    %         hpr = -6:0.5:15;
    %         obj = zeros(size(beta,2),size(rho,2),size(hpr,2));
    %         for nj = 1:size(beta,2)
    %             for nr = 1:size(rho,2)
    %                 for nm = 1:size(hpr,2)
    %                     hpstart = [beta(nj) rho(nr) hpr(nm)]';
    %                     obj(nj,nm) = nglglklhd_simp(hpstart,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)],lbi,ubi);
    %                 end
    %             end
    %         end
    %
    %         [~, ind_bbeta] = min(obj(:));
    %         [indb, indr, indc] = ind2sub([nj nr nm],ind_bbeta);
    %         hptmp = [beta(indb) rho(indr) hpr(indc)]';
    %     end
    
    
    if strcmp(kernel{ni},'NEW2')
        
        %     beta =[0.5:0.1:0.9 0.99]';
        %     hpr = [-6 -4 -2 0:1:5];
        beta =[0.5:0.01:0.9 0.95];
        hpr = -6:0.5:15;
        betahpr = beta;
        obj = zeros(size(beta,2),size(betahpr,2),size(hpr,2));
        for nj = 1:size(beta,2)
            for nr = 1:size(betahpr,2)
                for nm = 1:size(hpr,2)
                    hpstart = [beta(nj) betahpr(nr) hpr(nm)]';
                    obj(nj,nr,nm) = nglglklhd_simp(hpstart,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)],lbi,ubi);
                end
            end
        end
        [~, ind_bbeta] = min(obj(:));
        [indj, indr, indc] = ind2sub([nj nr nm],ind_bbeta);
        hptmp = [beta(indj) betahpr(indr) hpr(indc)]';
    end
    
    if strcmp(kernel{ni},'SSP')
        beta = 0.5:0.03:0.96;
        betahpr = 0.1:0.03:0.96;
        gamma = (0:0.05:1)*2*pi;
        hpr = -6:1:15;
        
        obj = zeros(size(beta,2),size(betahpr,2),size(gamma,2),size(hpr,2));
        for nj = 1:size(beta,2)
            for nr = 1:size(betahpr,2)
                for ng = 1:size(gamma,2)
                    for nm = 1:size(hpr,2)
                        hpstart = [beta(nj) betahpr(nr) gamma(ng) hpr(nm)]';
                        try
                            obj(nj,nr,ng,nm) = nglglklhd_simp(hpstart,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)],lbi,ubi);
                        catch
                            obj(nj,nr,ng,nm) = +inf;
                        end
                    end
                end
            end
        end
        [~, ind_bbeta] = min(obj(:));
        [indj, indr, indc, indd] = ind2sub([nj nr ng nm],ind_bbeta);
        hptmp = [beta(indj) betahpr(indr) gamma(indc) hpr(indd)]';
        
    end
    
    if strcmp(kernel{ni},'ML1')
        %lambda alpha(freq) c
        %     beta =[0.5:0.1:0.9 0.99]';
        %     hpr = [-6 -4 -2 0:1:5];
        beta =[0.5:0.01:0.9 0.95];
        betahpr = [log([1e-5 1e-4 1e-3 0.01:0.05:1]) log(1.1:0.1:2) 1 1.5 2];
        hpr = -6:0.5:15;

        obj = zeros(size(beta,2),size(betahpr,2),size(hpr,2));
        for nj = 1:size(beta,2)
            for nr = 1:size(betahpr,2)
                for nm = 1:size(hpr,2)
                    hpstart = [beta(nj) betahpr(nr) hpr(nm)]';
                    obj(nj,nr,nm) = nglglklhd_simp(hpstart,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)],lbi,ubi);
                end
            end
        end
        [~, ind_bbeta] = min(obj(:));
        [indj, indr, indc] = ind2sub([nj nr nm],ind_bbeta);
        hptmp = [beta(indj) betahpr(indr) hpr(indc)]';
    end
    
    if strcmp(kernel{ni},'ML2')
        
        beta =[0.5:0.01:0.9 0.95];
        betahpr = 2*pi*[0:0.05:1];
        hpr = -6:0.5:1;
        gamma = 0.99;
        
        obj = zeros(size(beta,2),size(betahpr,2),size(gamma,2),size(hpr,2));
        for nj = 1:size(beta,2)
            for nr = 1:size(betahpr,2)
                for ng = 1:size(gamma,2)
                    for nm = 1:size(hpr,2)
                        hpstart = [beta(nj) betahpr(nr) gamma(ng) hpr(nm)]';
                        try
                            obj(nj,nr,ng,nm) = nglglklhd_simp(hpstart,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)],lbi,ubi);
                        catch
                            obj(nj,nr,ng,nm) = +inf;
                        end
                    end
                end
            end
        end
        [~, ind_bbeta] = min(obj(:));
        [indj, indr, indc, indd] = ind2sub([nj nr ng nm],ind_bbeta);
        hptmp = [beta(indj) betahpr(indr) gamma(indc) hpr(indd)]';
        
    end

    if strcmp(kernel{ni},'ML3')
        
        beta =[0.5:0.01:0.9 0.95];
        betahpr = [log([1e-5 1e-4 1e-3 0.01:0.05:1]) log(1.1:0.1:2) 1 1.5 2];
        hpr = -6:0.5:15;
        gamma = [log([1e-5 1e-4 1e-3 0.01:0.05:1]) log(1.1:0.1:2) 1 1.5 2];
        
        obj = zeros(size(beta,2),size(betahpr,2),size(gamma,2),size(hpr,2));
        for nj = 1:size(beta,2)
            for nr = 1:size(betahpr,2)
                for ng = 1:size(gamma,2)
                    for nm = 1:size(hpr,2)
                        hpstart = [beta(nj) betahpr(nr) gamma(ng) hpr(nm)]';
                        try
                            obj(nj,nr,ng,nm) = nglglklhd_simp(hpstart,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)],lbi,ubi);
                        catch
                            obj(nj,nr,ng,nm) = +inf;
                        end
                    end
                end
            end
        end
        [~, ind_bbeta] = min(obj(:));
        [indj, indr, indc, indd] = ind2sub([nj nr ng nm],ind_bbeta);
        hptmp = [beta(indj) betahpr(indr) gamma(indc) hpr(indd)]';
        
    end
    
    if strcmp(kernel{ni},'SI2')
        beta = 0.5:0.03:0.96;
        betahpr = 0.1:0.05:0.95;
        gamma = [0.1:0.1:0.9 0.95];
        hpr = -6:1:1;
        
        obj = zeros(size(beta,2),size(betahpr,2),size(gamma,2),size(hpr,2));
        for nj = 1:size(beta,2)
            for nr = 1:size(betahpr,2)
                for ng = 1:size(gamma,2)
                    for nm = 1:size(hpr,2)
                        hpstart = [beta(nj) betahpr(nr) gamma(ng) hpr(nm)]';
                        try
                            obj(nj,nr,ng,nm) = nglglklhd_simp(hpstart,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)],lbi,ubi);
                        catch
                            obj(nj,nr,ng,nm) = +inf;
                        end
                    end
                end
            end
        end
        [~, ind_bbeta] = min(obj(:));
        [indj, indr, indc, indd] = ind2sub([nj nr ng nm],ind_bbeta);
        hptmp = [beta(indj) betahpr(indr) gamma(indc) hpr(indd)]';
        
    end
    
    if Nu > 1
        hpstart = hptmp;
        ff = @(x)nglglklhd(x,Rti,Ney,kernel(ni),sigma,[0 nb(ni)+1-nk(ni)], [0 length(hpstart)]);
        hptmp = fmincon(ff,hpstart,[],[],[],[],lbi,ubi,[],options);
    end
    
    hpini = [hpini; hptmp];
    lb = [lb; lbi];
    ub = [ub; ubi];
    ind2(ni) = length(hpstart);
    
end
ind2 = cumsum(ind2); ind2=[0 ind2];



function [Rt, Ney, sigma] = qrfactor_Phi(z,nb,nk)

[Nd, ~, Nu, Ne] = size(z);
ind = cumsum(nb+1-nk); ind = [0 ind]; % index of the estimated FIR parameters for each of Nu blocks
Nth = ind(end);
Rt = zeros(0,Nth+1);


dom = pvget(z,'Domain');

if lower(dom(1)) == 't'
    
    y = pvget(z,'OutputData');
    u = pvget(z,'InputData');
    
    Ney = sum(Nd-max(nb));
    yt = zeros(Ney,1);
    Phi = zeros(Nth,Ney);
    nt = cumsum(Nd-max(nb)); nt = [0 nt];
    
    for ni = 1:Nu
        for ne = 1:Ne
            for nj = 1:Nd(ne)-max(nb)
                Phi(ind(ni)+1:ind(ni+1),nj+nt(ne)) = flipud(u{ne}(nj+max(nb)-nb(ni):(nj+max(nb)-nk(ni)),ni));
            end
        end
    end
    
    
    for ne = 1:Ne
        yt(nt(ne)+1:nt(ne+1)) = y{ne}(max(nb)+1:end);
    end
    
    % Calculate the QR factor of Phi and yt
    Rt = triu(qr([Phi' yt]));
    try
        Rt = Rt(1:Nth+1,:);
    catch
        Rt = [Phi' yt];
    end
    
end

if lower(dom(1)) == 'f'
    
    if min(z.Frequency) == 0
        zt = ifft(z);
        z = fft(zt,size(zt,1),'compl');
    end
    
    nz = nyqcut(complex(z));
    y = pvget(nz,'OutputData');
    u = pvget(nz,'InputData');
    w = pvget(nz,'Radfreqs');
    Ts = pvget(nz,'Ts');
    
    Nd = size(nz,1);
    Ney = sum(Nd);
    
    nt = cumsum(Nd); nt = [0 nt];
    yt = zeros(Ney,1);
    Phi = zeros(Nth,Ney);
    
    for ne = 1:Ne
        for ni = 1:Nu
            OM = exp(-1i*(nk(ni):nb(ni))'*w{ne}'*Ts{ne});
            tempM = (OM.').*(u{ne}(:,ni)*ones(1,nb(ni)-nk(ni)+1));
            Phi(ind(ni)+1:ind(ni+1),nt(ne)+1:nt(ne+1)) = tempM.';
        end
        yt(nt(ne)+1:nt(ne+1)) = y{ne};
    end
    
    % Calculate the QR factor of Phi and yt
    if norm(imag(triu(qr([Phi.' yt])))) > 1e-4
        warning('The frequency domain data processing is not handled appropriately!')
    end
    Rt = real(triu(qr([Phi.' yt])));
    try
        Rt = Rt(1:Nth+1,:);
    catch
        Rt = [Phi.' yt];
    end
    
end

% Calculate the variance of the measurement noise

if Ney > Nth && cond((Rt(:,1:Nth)'*Rt(:,1:Nth)))<1e+12
    X = (Rt(:,1:Nth)'*Rt(:,1:Nth))\Rt(:,1:Nth)'*Rt(:,end);
    sigma = sum((Rt(:,end)-Rt(:,1:Nth)*X).^2)/(Ney-Nth);
else
    Nrdth = floor(min(0.8*Ney,Nth)/Nu);
    rdthind = [];
    for ni = 1:Nu
        rdthind = [rdthind ind(ni)+1:ind(ni)+min(nb(ni),Nrdth)];
    end
    Nrdy = 0;
    X = (Rt(Nrdy+1:end,rdthind)'*Rt(Nrdy+1:end,rdthind))\Rt(Nrdy+1:end,rdthind)'*Rt(Nrdy+1:end,end);
    sigma = sum((Rt(Nrdy+1:end,end)-Rt(Nrdy+1:end,rdthind)*X).^2)/(Ney-Nrdy-length(rdthind));
end






function [obj, condNb, xest] = nglglklhd_simp(hp,Rt,Ney, kernel,sigma,ind,ind2,lb,ub)

%warning('off', 'MATLAB:nearlySingularMatrix');
warning off;

Nth = size(Rt,2) -1;
Nu = size(kernel,2);
P_chol = zeros(Nth,Nth);


chol_flag = 1;
cond_flag = 1;
ini_flag = 1;


for di = 1:Nu
    
    hyper = hp(ind2(di)+1:ind2(di+1));
    nPi = ind(di+1)-ind(di);
    Pi = zeros(nPi,nPi);
    
    lbi = lb(ind2(di)+1:ind2(di+1));
    ubi = ub(ind2(di)+1:ind2(di+1));
    
    ini_flag_di = sum(hyper<=ubi) + sum(hyper>=lbi);
    if ini_flag_di ~= 2*length(hyper)
        ini_flag = 0;
        break;
    end
    
    switch kernel{di}
        
        case 'CS'
            for dk=1:nPi
                for dj=1:nPi
                    Pi(dk,dj) = dk*dj*min(dk,dj)/2-min(dk,dj)^3/6;
                end
            end
            
        case 'SE'
            for dk=1:nPi
                for dj=1:nPi
                    Pi(dk,dj) = exp(-(dk-dj)^2/(2*hyper(1)^2));
                end
            end
            
            
        case 'DI'
            t=abs(hyper(1)).^(1:nPi);%exp(-beta*[1:nb])
            Pi=diag(t);
            
        case 'SS'
            t=abs(hyper(1)).^(1:nPi);%exp(-beta*[1:nb])
            for dk=1:nPi
                for dj=1:nPi
                    Pi(dk,dj) = t(dk)*t(dj)*min(t(dk),t(dj))/2-min(t(dk),t(dj))^3/6;
                end
            end
        case 'SSP'
            t=abs(hyper(1)).^(1:nPi);%exp(-beta*[1:nb])
            for dk=1:nPi
                for dj=1:dk
                    Pi(dk,dj) = t(dk)*t(dj)*min(t(dk),t(dj))/2-min(t(dk),t(dj))^3/6;
                    Pi(dj,dk) = Pi(dk,dj);
                end
            end            
            a = dimpulse([1 zeros(1,2)],[1 -2*hyper(2)*cos(hyper(3)), hyper(2)^2],nPi);
            B = toeplitz(a,zeros(1,nPi));
            Pi = B*Pi*B';
            
        case 'HF'
            t=abs(hyper(1)).^(1:nPi);%exp(-beta*[1:nb])
            for dk=1:nPi
                for dj=1:nPi
                    if mod(dk+dj,2) == 0
                        Pi(dk,dj) = min(t(dk),t(dj));
                    else
                        Pi(dk,dj) = -min(t(dk),t(dj));
                    end
                end
            end
            
        case 'TC'
            t=abs(hyper(1)).^(1:nPi); %exp(-beta*[1:nb])
            for dk=1:nPi
                for dj=1:dk
                    Pi(dk,dj) = min(t(dk),t(dj));
                    Pi(dj,dk) = Pi(dk,dj);
                end
            end
            
        case 'DC'
            for dk = 1:nPi
                for dj = 1:dk
                    Pi(dk,dj) = abs(hyper(1))^((dk+dj)/2)*hyper(2)^(abs(dk-dj));
                    Pi(dj,dk) = Pi(dk,dj); 
                end
            end
        case 'NEW1'
            %             for dk = 1:nPi
            %                 for dj = 1:nPi
            %                     %                 Pi(dk,dj) = min(abs(hyper(1))^dk,abs(hyper(1))^dj)-abs(hyper(1))^(dk+dj);
            %                     Pi(dk,dj) = min(dk,dj)*abs(hyper(1))^(dk+dj);
            %                 end
            %             end
            for dk = 1:nPi
                for dj = 1:nPi
                    Pi(dk,dj) = (hyper(1))^max(dk,dj)+(hyper(1))^(dk+dj);
                end
            end
            %             for dk = 1:nPi
            %                 for dj = 1:nPi
            %                     Pi(dk,dj) = (1-(min(dk,dj)+1)^(-1))*(hyper(1))^(dk+dj);
            %                 end
            %             end
            
        case 'NEW11'
            %             for dk = 1:nPi
            %                 for dj = 1:nPi
            %                     %                 Pi(dk,dj) = min(abs(hyper(1))^dk,abs(hyper(1))^dj)-abs(hyper(1))^(dk+dj);
            %                     Pi(dk,dj) = min(dk,dj)*abs(hyper(1))^(dk+dj);
            %                 end
            %             end
            %             for dk = 1:nPi
            %                 for dj = 1:nPi
            %                     Pi(dk,dj) = (hyper(1))^max(dk,dj)+(hyper(1))^(dk+dj);
            %                 end
            %             end
            for dk = 1:nPi
                for dj = 1:nPi
                    Pi(dk,dj) = (1-(min(dk,dj)+1)^(-1))*(hyper(1))^(dk+dj);
                end
            end
            
        case 'NEW2'
            for dk = 1:nPi
                for dj = 1:nPi
                    if hyper(1) == hyper(2)
                        %                 Pi(dk,dj) = min(abs(hyper(1))^dk,abs(hyper(1))^dj)-abs(hyper(1))^(dk+dj);
                        Pi(dk,dj) = min(dk,dj)*(hyper(1))^(dk+dj);
                    else
                        Pi(dk,dj) = sign(hyper(2)-hyper(1))*((hyper(1))^(max(dk,dj)-min(dk,dj))*(hyper(2))^(2*min(dk,dj)) - (hyper(1))^(dk+dj));
                    end
                end
            end
            
        case 'ML1'
            for dk = 1:nPi
                for dj = 1:dk
                        Pi(dk,dj) = hyper(1)^(dk+dj)*cos(exp(hyper(2))*abs(dk-dj));
                        Pi(dj,dk) = Pi(dk,dj);
                end
            end
            
        case 'ML2'
            for dk = 1:nPi
                for dj = 1:dk
                        Pi(dk,dj) = hyper(1)^(dk+dj)*(sin(hyper(2)*dk)+1)*(sin(hyper(2)*dj)+1)*hyper(3)^(abs(dk-dj));
                        Pi(dj,dk) = Pi(dk,dj);
                end
            end
        case 'ML3'
            for dk = 1:nPi
                for dj = 1:dk
                        Pi(dk,dj) = hyper(1)^(dk+dj)*(sin(exp(hyper(2))*dk)+1)*(sin(exp(hyper(2))*dj)+1)*cos(exp(hyper(3))*abs(dk-dj));
                        Pi(dj,dk) = Pi(dk,dj);
                end            
            end
        case 'SI2'
            %             for dk = 1:nPi
            %                 for dj = 1:nPi
            %                     %                     Pi(dk,dj) = (hyper(1))^(max(dk,dj)) + (exp(hyper(2))-1)*(hyper(1))^(dk+dj);
            %                     Pi(dk,dj) = (1+exp(hyper(2))-(min(dk,dj)+1)^(-1))*(hyper(1))^(dk+dj);
            %                 end
            %             end
            % hyper = [e^(-w0*xi), xi(damping ratio), e^(-gamma) uncertainty, scaling rate]
            be = -log(hyper(1))/hyper(2)*sqrt(1-hyper(2)^2);
            theta = acos(2*log(hyper(3)/hyper(1))/sqrt(4*be^2+4*log(hyper(3)/hyper(1))^2));
            for dk = 1:nPi
                for dj = 1:dk
                    Pi(dk,dj) = hyper(1)^(dk+dj)*(cos(be*dk)-log(hyper(1))/be*sin(be*dk))*(cos(be*dj)-log(hyper(1))/be*sin(be*dj)) ...
                        + hyper(1)^(dk+dj)*sin(be*dk)*sin(be*dj)/be^2 + hyper(1)^(dk+dj)*cos(be*(dk-dj))*((hyper(3)/hyper(1))^(2*min(dk,dj))-1)/(4*be^2*log(hyper(3)/hyper(2)))...
                        + hyper(1)^(dk+dj)*(cos(theta+be*(dk+dj)) - (hyper(3)/hyper(1))^(2*min(dk,dj))*cos(2*be*min(dk,dj)-theta-be*(dk+dj)) )...
                        /(2*be^2*sqrt(4*be^2+4*log(hyper(3)/hyper(2))^2));
                    Pi(dj,dk) = Pi(dk,dj);
                end
            end
    end
    
    if sum(sum(isinf(Pi))) ~= 0 || sum(sum(isnan(Pi))) ~= 0
        cond_flag = 0;
        break;
    end
    
    if cond(Pi) > 1e+100
        cond_flag = 0;
    end
    
    
    try
        U = chol(Pi)'*sqrt(exp(hyper(end)));
    catch
        try
            U=chol(Pi+1e-4*eye(nPi))'*sqrt(exp(hyper(end)));
        catch
            chol_flag = 0;
            break;
        end
    end
    P_chol(ind(di)+1:ind(di+1),ind(di)+1:ind(di+1)) = U;
    
end

if ini_flag ~= 1 || cond_flag ~= 1 || chol_flag ~=1
    obj = inf;
    return;
end

R1 = real(triu(qr([Rt(:,1:Nth)*P_chol Rt(:,Nth+1); sqrt(sigma)*eye(Nth) zeros(Nth,1)])));
R1 = R1(1:Nth+1,:);

obj = R1(end,Nth+1)^2/sigma + 2*sum(log(abs(diag(R1(1:Nth,1:Nth)))))+ (Ney-Nth)*log(sigma);

if nargout >1
    condNb = cond((R1(1:Nth,1:Nth)'*R1(1:Nth,1:Nth)));
    xest = P_chol/R1(1:Nth,1:Nth)*R1(1:Nth,Nth+1);
end


