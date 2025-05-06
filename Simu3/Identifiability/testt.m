value = -2:0.03:8;

COST = zeros(length(value),1);
EFIT = zeros(length(value),1);
for ii = 1:length(value)

%  hpp = [value(ii);hp(2:end)];
% hpp = [hp(1);value(ii);hp(3:end)];
%  hpp = [hp(1:2);value(ii);hp(4:end)];
% hpp = [hp(1:3);value(ii);hp(5:end)];
% hpp = [hp(1:4);value(ii);hp(6:end)];
hpp = [hp(1:5);value(ii)];


[cost, Mx1, Mx2] = nglglklhd(hpp, Psi, y, kernel, M, method);
COST(ii) = cost;

% h0 = hp(end);
h0 = 0;
poly = hpp(1:M);

B = zeros(N-n+1,M);
for i = 1:M
    B(:,i) = u(n:N).^i;
end
b = B*poly;

switch method
    case 'chol'
        Oiinv =  Mx1'*Mx1;
    case 'svd'
        Oiinv =  Mx1*Mx2*Mx1';
    case 'combined'
        if Mx2(2,1) == 0
            Oiinv =  Mx1*Mx2*Mx1';
        else
            Oiinv =  Mx1'*Mx1;
        end
end


W = Oiinv*(y-h0-b);
[O, K, O1h] = CalculateOutputKernel(Psi, Psi, M, kernel, hpp, 1);
yhat = O*W + h0 + b;

d = load(['Databank/data_N' int2str(1000) '_repi=' int2str(37) '.mat']);
data = d.datainfo.data(1:N,:);
ytrue = d.datainfo.ytrue(n:N);
efit= gof(ytrue,yhat);
EFIT(ii) = efit;
end

figure('position', [500, 200, 1200, 600]); 
subplot(1,3,1);
plot(value,EFIT); grid on;
xlabel('\sigma^2'); ylabel('Efit');
subplot(1,3,2);
plot(value,COST); grid on;
xlabel('\sigma^2'); ylabel('cost');
subplot(1,3,3);
plot(COST,EFIT); grid on;
xlabel('cost'); ylabel('Efit');

% save('hp_ratio10.mat','hp')


% noise = randn(1000,1);
% noise = (noise - mean(noise))*std(datainfo.ytrue)/std(noise)/sqrt(datainfo.snr);
% datainfo.noise = noise;
% datainfo.NoiseVariance = var(noise);
% datainfo.data(:,2) = datainfo.ytrue + datainfo.noise;
% 
% save('data_N1000_repi=37_snr10.mat','datainfo')