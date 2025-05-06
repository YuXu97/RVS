
subplot(3,1,1)
PFIT_M2 = [PFIT_M2s1 PFIT_M2s5 PFIT_M2s10 PFIT_M2s20]; 
nice_groupboxplot(PFIT_M2,20,100,{'SED-MPK','DC-decay-w','DC-ob-w',},{'SNR=1dB','SNR=5dB','SNR=10dB','SNR=20dB'}, ...
    2, 'M=2','',2)

subplot(3,1,2)
PFIT_M3 = [PFIT_M3s1 PFIT_M3s5 PFIT_M3s10 PFIT_M3s20]; 
nice_groupboxplot_nolg(PFIT_M3,20,100,{'SED-MPK','DC-decay-w','DC-ob-w',},{'SNR=1dB','SNR=5dB','SNR=10dB','SNR=20dB'}, ...
    2, 'M=3','',2)

subplot(3,1,3)
PFIT_M4 = [PFIT_M4s1 PFIT_M4s5 PFIT_M4s10 PFIT_M4s20]; 
nice_groupboxplot_nolg(PFIT_M4,20,100,{'SED-MPK','DC-decay-w','DC-ob-w',},{'SNR=1dB','SNR=5dB','SNR=10dB','SNR=20dB'}, ...
    2, 'M=4','',2)