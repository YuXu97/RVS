clear;
path = [pwd '/exp3_1/M3n80N400L30SNR10dB_Gaussianinput/'];
snr = db2pow(5);
Maxrepi = 200;
Nmax = 4000;

for repi = 1:Maxrepi
    d = load([path '/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat']);
    datainfo.ytrue = d.datainfo.ytrue;
    datainfo.LinearSystem1 = d.datainfo.LinearSystem1;
    datainfo.LinearSystemPoles1 = d.datainfo.LinearSystemPoles1;
    datainfo.LinearSystemImpulseResponse1 = d.datainfo.LinearSystemImpulseResponse1;
    datainfo.LinearSystem2 = d.datainfo.LinearSystem2;
    datainfo.LinearSystemPoles2 = d.datainfo.LinearSystemPoles2;
    datainfo.LinearSystemImpulseResponse2 = d.datainfo.LinearSystemImpulseResponse2;
    datainfo.h0 = d.datainfo.h0;
    datainfo.LinearSysOrder = d.datainfo.LinearSysOrder;
    datainfo.PolyCoefficients = d.datainfo.PolyCoefficients;
    %datainfo.HmVar = d.datainfo.HmVar;
    datainfo.snr = snr;
    
    noise = randn(Nmax,1);
    ytrue = d.datainfo.ytrue;
    noise = (noise - mean(noise))*std(ytrue)/std(noise)/sqrt(snr);
    data = d.datainfo.data;
    data = [data(:,1) ytrue+noise];
    datainfo.data = data;
    
    datainfo.noise = noise;
    datainfo.NoiseVariance = var(noise);
    
    return_filename = [pwd '/data_N' int2str(Nmax) '_repi=' int2str(repi) '.mat'];
    save(return_filename, 'datainfo');
end


