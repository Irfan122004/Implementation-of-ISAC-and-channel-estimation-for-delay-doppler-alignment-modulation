
clear; clc; close all;

%% Transmit power range
P_dBm = 20:5:50;
P_noise_dBm = -94;

% Convert to Watts
P_noise = 10^(P_noise_dBm/10)/1000;
P_t = 10.^(P_dBm/10)/1000;

%% System parameters
L = 5;              % number of paths
N = 10000;          % total symbols in coherence block
N_p = 100;          
N_g = 100;
N_d = N - N_p - 2*N_g;   % data symbols

%% Channel (fixed realization for fairness)
h_gain = abs(randn(L,1) + 1j*randn(L,1)).^2;
f_gain = ones(L,1);   % ZF beamforming gain

%% DDAM-ZF SNR
SNR_comm = (sum(sqrt(h_gain .* f_gain)).^2) .* P_t / P_noise;

%% ✅ Correct DDAM Spectral Efficiency
R_DDAM_ZF = (N_d / N) * log2(1 + SNR_comm);

%% OFDM parameters
W = 512;            % subcarriers
N_cp = 100;         % cyclic prefix
N_OFDM_sym = floor(N_d / (W + N_cp));

%% OFDM SNR per subcarrier
SNR_per_sub = SNR_comm / W;

%% OFDM Spectral Efficiency
R_OFDM = (N_OFDM_sym * W / N) * log2(1 + SNR_per_sub);

%% Plot
figure;
plot(P_dBm, R_DDAM_ZF, 'b-o', 'LineWidth', 2); hold on;
plot(P_dBm, R_OFDM, 'r-s', 'LineWidth', 2);

xlabel('Transmit Power (dBm)');
ylabel('Spectral Efficiency (bits/s/Hz)');
legend('DDAM-ZF', 'OFDM', 'Location', 'NorthWest');
title('Spectral Efficiency Comparison (Corrected)');
grid on;

%% Improve visualization
ylim([0 max(R_DDAM_ZF)*1.2]);