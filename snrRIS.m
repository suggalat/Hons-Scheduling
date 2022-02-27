function [SNR]=snrRIS(V_ch,w)
    load('ga.mat');
    SNR = 10*log10((P/B)*abs(V_ch*w)^2/noise_pwr);
    end