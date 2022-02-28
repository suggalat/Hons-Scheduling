function [fitVal]=fitness(w)
    load('ga.mat');
    snr1=snrRIS(Ch_Rx1,w);
    snr2=snrRIS(Ch_Rx2,w);
%     fitVal= snr1+snr2+200; % Offset - 200
    fitVal=1/pdist([snr1 snr2;20 20]);
end