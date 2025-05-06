data_M2 = [PFIT_M2snr1 PFIT_M2snr5 PFIT_M2snr10];
data_M3 = [PFIT_M3snr1 PFIT_M3snr5 PFIT_M3snr10];

subplot(2,1,1)
nice_groupboxplot(data_M2,20,100,{'SED-MPK','DC-decay-w','DC-ob-w'},{'SNR=1dB','SNR=5dB','SNR=10dB'}, 1, 'M=2','',2);

subplot(2,1,2)
nice_groupboxplot_nolg(data_M3,20,100,{'SED-MPK','DC-decay-w','DC-ob-w'},{'SNR=1dB','SNR=5dB','SNR=10dB'}, 1, 'M=3','',2);

p_boxplot(data_M3,0,100,{'SED-MPK','DC-decay-w','DC-ob-w','SED-MPK','DC-decay-w','DC-ob-w','SED-MPK','DC-decay-w','DC-ob-w'},'','')
