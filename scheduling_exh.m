% Execute this program, only after executing scheduling.m
% This program inputs G,H and other parameters obtained in scheduling.m
Tx_xyz = [1,1,1];           % Tx co-ordinates
Rx1_xyz = [30,1,1];         % Rx-1 co-ordinates
Rx2_xyz = [29,3,1];         % Rx-2 co-ordinates
RIS_xyz = [25,6,2];         % RIS co-ordinates
K = 2;                      % Number of Subcarriers

SNR = 5;                    % SNR Without RIS in dB
n_fft = 2;
n_cpe = 0;
Frequency = 6;              % Frequency in GHz
ArrayType = 2;              % Uniform Linear Array=1 or Uniform Planar Array=2
Environment = 1;            % 1 Indoor (InH - Indoor Office) / 2 Outdoor (UMi - Street Canyon)
% N = 16;                      % Number of RIS elements
Nsym = 1000;                % Number of Realisations
Nt = 1;                     % Number of antennas at Tx
Nr = 1;                     % Number of antennas at Rx
Scenario=1;                 % 1 (RIS in xz plane - left side wall) or 2 (RIS in yz plane - opposite wall)
% Had = hadamard(N);
% Had_flip = Had;
% Had_flip(1,:) = -1;
% RIS_config = [Had, Had_flip, -Had, -Had_flip];              % 4N possible RIS configurations

C = cell(1,N);
[C{:}] = ndgrid([1,-1]);
C = cellfun(@(a)a(:),C,'Uni',0);
M = [C{:}];

RIS_config = M';   % 'N' RIS configurations used for channel estimation
H1 = zeros(N,Nt,Nsym);
G1 = zeros(Nr,N,Nsym);
D  = zeros(Nr,Nt,Nsym);
H2 = zeros(N,Nt,Nsym);
P = 1;  % Power in Watts
B = 10e6; % bandwidth
mod_method = 'QPSK';

Noise_iter=100;
% Calculate modulation order from modulation method
mod_methods = {'BPSK','QPSK','8PSK','16QAM','32QAM','64QAM'};
mod_order = find(ismember(mod_methods,mod_method));

%% Symbol modulation
% BPSK
if mod_order == 1
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n);
    quadrature = sin(n);
    symbol_book = (in_phase + quadrature*1i)';
end

% Phase shift keying about unit circle
if mod_order == 2 || mod_order == 3
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n+pi/4);
    quadrature = sin(n+pi/4);
    symbol_book = (in_phase + quadrature*1i)';
end

% 16QAM, 64QAM modulation 
if mod_order == 4 || mod_order == 6
   mod_ind = sqrt(2^mod_order);
   in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
   quadrature = repmat(linspace(-1,1,mod_ind)',1,mod_ind);
   symbol_book = in_phase(:) + quadrature(:)*1i;
end

% 32QAM modulation
% Generates 6x6 constallation and removes corners
if mod_order == 5
   mod_ind = 6;
   in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
   quadrature = repmat(linspace(-1,1,mod_ind)',1,mod_ind);
   symbol_book = in_phase(:) + quadrature(:)*1i;
   symbol_book = symbol_book([2:5 7:30 32:35]);
end

X =sqrt(P/B)* [symbol_book(1);symbol_book(3)]; % Modulate data according to symbol_book
X = repmat(X,1,length(RIS_config));
            
%% Use IFFT to move to time domain
% Pad input signal to appropriate length
fft_rem = mod(n_fft-mod(length(X),n_fft),n_fft);
X_padded = [X;zeros(fft_rem,1)];
X_blocks = reshape(X_padded,n_fft,[]);
x = ifft(X_blocks);

% Add cyclic prefix extension and shift from parellel to serial
x_cpe = [x(end-n_cpe+1:end,:);x];
x_s = x_cpe(:);

%% SimRIS

% for iter = 1:Nsym  
%    [H1(:,:,iter),G11,D(:,:,iter)] = SimRIS_v18_1(Environment,Scenario,Frequency,ArrayType,N,Nt,Nr,Tx_xyz,Rx1_xyz,RIS_xyz);
%    [H2(:,:,iter),G12,D(:,:,iter)] = SimRIS_v18_1(Environment,Scenario,Frequency,ArrayType,N,Nt,Nr,Tx_xyz,Rx2_xyz,RIS_xyz);   
%    G1(:,:,iter) = (G11+G12)/2;
% end
% G1 = mean(G1,3);
% H1 = mean(H1,3);
% H2 = mean(H2,3);
% Ch_Rx1 = (G1'.*H1)';
% Ch_Rx2 = (G1'.*H2)';

%% Add AWGN

x_s_ch1 = (Ch_Rx1*RIS_config).*x_cpe;
x_s_ch2 = (Ch_Rx2*RIS_config).*x_cpe;

% Calculate data power
data_pwr = mean(abs(x_s_ch1(:,1).^2));

%Add noise to channel
noise_pwr = data_pwr/10^(SNR/10);

[sizeX_row,sizeX_col]=size(x_s_ch1);
noise = zeros(sizeX_row,sizeX_col,Noise_iter);
x_s_noise1=zeros(sizeX_row,sizeX_col,Noise_iter);
x_s_noise2=zeros(sizeX_row,sizeX_col,Noise_iter);
SNR1=zeros(sizeX_col,Noise_iter);
SNR2=zeros(sizeX_col,Noise_iter);
for iter = 1:Noise_iter
    noise(:,:,iter) = normrnd(0,sqrt(noise_pwr/2),size(x_s_ch1)) + normrnd(0,sqrt(noise_pwr/2),size(x_s_ch1))*1i;
    x_s_noise1(:,:,iter)=x_s_ch1+noise(:,:,iter);
    SNR1(:,iter)=10*log10(abs(mean(x_s_ch1.^2)./mean(noise(:,:,iter).^2)));
    x_s_noise2(:,:,iter)=x_s_ch2+noise(:,:,iter);
    SNR2(:,iter)=10*log10(abs(mean(x_s_ch2.^2)./mean(noise(:,:,iter).^2)));
end

SNR_Rx1 = mean(SNR1,2);
SNR_Rx2 = mean(SNR2,2);
x_s_noise1=mean(x_s_noise1,3);
x_s_noise2=mean(x_s_noise2,3);

% x_s_noise1 = x_s_ch1 + mean(noise,3);
% x_s_noise2 = x_s_ch2 + mean(noise,3);

% snr_mean1 = 10*log10(abs(mean(x_s_ch1.^2)/mean(noise.^2)));
% % snr_mean1 = 10*log10(mean(abs(x_s_ch1.^2))/mean(abs(noise1.^2)));
% % snr_mean2 = 10*log10(mean(abs(x_s_ch2.^2))/mean(abs(noise2.^2)));
% snr_mean2 = 10*log10(abs(mean(x_s_ch2.^2)/mean(noise.^2)));

%% Use FFT to move to frequency domain
% Remove cyclic prefix extension and shift from serial to parellel
x_p1 = reshape(x_s_noise1,n_fft+n_cpe,[]);
x_p_cpr1 = x_p1(n_cpe+1:end,:);

% Move to frequency domain
X_hat_blocks1 = fft(x_p_cpr1);

% K1 = lsqr(X_blocks,X_hat_blocks1);
K1 = X_hat_blocks1./X_blocks;
Ch_Rx1_est = K1/RIS_config; 
Rx1_rec = X_hat_blocks1./K1; % Received Rx1
err1 = mean(abs(Ch_Rx1_est-Ch_Rx1));

x_p2 = reshape(x_s_noise2,n_fft+n_cpe,[]);
x_p_cpr2 = x_p2(n_cpe+1:end,:);

% Move to frequency domain
X_hat_blocks2 = fft(x_p_cpr2);
% K2 = lsqr(X_blocks,X_hat_blocks2);
K2 = X_hat_blocks2./X_blocks;
Ch_Rx2_est = K2/RIS_config;
Rx2_rec = X_hat_blocks2./K2;
err2 = mean(abs(Ch_Rx2_est-Ch_Rx2));

% SNR_Rx1 = 10*log10(abs(mean(x_p1.^2))/noise_pwr);
% SNR_Rx2 = 10*log10(abs(mean(x_p2.^2)/noise_pwr));
% SNR_Rx1 = (10*log10((P/(B*noise_pwr))*(mean(abs(K1)).^2)));
% SNR_Rx2 = (10*log10((P/(B*noise_pwr))*(mean(abs(K2)).^2)));
max_SNR = max(max(SNR_Rx1),max(SNR_Rx2)); 
min_SNR=min(min(SNR_Rx1),min(SNR_Rx2));


temp = [SNR_Rx1,SNR_Rx2];
positive_coords=temp(all((temp-5)>0,2),:);
[farVal,farInd]=max(pdist2([0,0],positive_coords));
selectedPoint = positive_coords(farInd,:);
if(length(positive_coords)==0)
    dist_Th=pdist2([SNR,SNR],[SNR_Rx1,SNR_Rx2]);
    [minVal,minIndex] = min(dist_Th);
    selectedPoint=[SNR_Rx1(minIndex),SNR_Rx2(minIndex)];
end
tempIn = find(selectedPoint==temp)
selectedConfig = RIS_config(:,tempIn(1));

% quiver(0,0,real(symbol_book(3)),complex(symbol_book(3)));
for iter=1:N
   Cons(iter)=Ch_Rx2(iter)*selectedConfig(iter)*symbol_book(3);
%    quiver(0,0,real(Cons(iter)),complex(Cons(iter)));   
end
symb=symbol_book(3);

scatter(real(symb),complex(symb));
hold on;
scatter(real(Cons),complex(Cons));
quiver(0,0,real(sum(Cons)),complex(sum(Cons)));
%%  PLOTS
figure
plot(1:N,(Ch_Rx1),1:N,(mean(Ch_Rx1_est)));grid on
legend('original','estimated')
xlabel('Index');
ylabel('Channel Coeff');
title(sprintf('\\bfChannel Estimation-Rx1\n\\rm Estimation Error %e',mean(err1)));

figure 
scatter(SNR_Rx1,SNR_Rx2,'filled');
xlabel('SNR achieved at Rx1');
ylabel('SNR achieved at Rx2');
hold on
plot(floor(min_SNR):ceil(max_SNR),SNR*ones(1,ceil(max_SNR)-floor(min_SNR)+1));
plot(SNR*ones(1,ceil(max_SNR)-floor(min_SNR)+1),floor(min_SNR):ceil(max_SNR));
scatter(selectedPoint(1),selectedPoint(2),'red','filled');
title(sprintf('\\bfSNR achieved at Rx2 vs Rx1'));

figure
plot(1:N,(Ch_Rx2),1:N,(mean(Ch_Rx2_est)));grid on
legend('original','estimated')
xlabel('Index');
ylabel('Channel Coeff');
title(sprintf('\\bfChannel Estimation-Rx2\n\\rm Estimation Error %e',mean(err2)));

figure
plot(1:length(RIS_config),SNR_Rx1,1:length(RIS_config),SNR_Rx2);
legend('SNR of Rx1','SNR of Rx2')
xlabel('Index');
ylabel('SNR(dB)');
title(sprintf('\\bf SNR of Rx1,Rx2 in dB '));


% C1 = X_hat_blocks1(:,1)./X_blocks(:,1);
% X_hat_blocks1 = X_hat_blocks1./repmat(C1,1,size(X_hat_blocks1,2));
% 
% C2 = X_hat_blocks2(:,1)./X_blocks(:,1);
% X_hat_blocks2 = X_hat_blocks2./repmat(C2,1,size(X_hat_blocks2,2));
