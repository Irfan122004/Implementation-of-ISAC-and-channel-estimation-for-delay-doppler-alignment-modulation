
clear; clc; close all;

K = 2048;           % OFDM subcarriers
L_vals = [5, 10, 30];
N_symbols = 10000;

% OFDM PAPR
ofdm_symbols = (randi(2, K, N_symbols)*2 - 1) + 1j*(randi(2, K, N_symbols)*2 - 1);
ofdm_papr = zeros(1, N_symbols);
for n = 1:N_symbols
    x = ifft(ofdm_symbols(:, n));
    ofdm_papr(n) = max(abs(x).^2) / mean(abs(x).^2);
end

% DDAM PAPR
ddam_papr = zeros(length(L_vals), N_symbols);
for l_idx = 1:length(L_vals)
    L = L_vals(l_idx);
    for n = 1:N_symbols
        s = (randi(2, L, 1)*2 - 1) + 1j*(randi(2, L, 1)*2 - 1);
        x = sum(s);
        ddam_papr(l_idx, n) = abs(x)^2 / mean(abs(s).^2);
    end
end

% CCDF
[ofdm_prob, ofdm_papr_db] = ecdf(10*log10(ofdm_papr));
ccdf_ofdm = 1 - ofdm_prob;
figure;
semilogy(ofdm_papr_db, ccdf_ofdm, 'b-', 'LineWidth', 1.5); hold on;
colors = {'r-', 'g-', 'm-'};
for l_idx = 1:length(L_vals)
    [prob, papr_db] = ecdf(10*log10(ddam_papr(l_idx, :)));
    semilogy(papr_db, 1-prob, colors{l_idx}, 'LineWidth', 1.5);
end
xlabel('PAPR (dB)'); ylabel('CCDF');
legend('OFDM (K=2048)', 'DDAM L=5', 'DDAM L=10', 'DDAM L=30');
grid on; title('PAPR Comparison');