clear; clc; close all;

J_vals = [5, 10, 20];
No_vals = [10, 100, 1000];
SNR_dB = -10:5:25;
N_iter = 500;
Tc = 0.1e-3;
nu_true = 500;

error = zeros(length(SNR_dB), length(J_vals), length(No_vals));

for idx_snr = 1:length(SNR_dB)
    SNR = 10^(SNR_dB(idx_snr)/10);
    noise_var = 1/SNR;
    for j_idx = 1:length(J_vals)
        J = J_vals(j_idx);
        for no_idx = 1:length(No_vals)
            No = No_vals(no_idx);
            err_sum = 0;
            for iter = 1:N_iter
                alpha = randn + 1j*randn;
                n = 0:J-1;
                phase = exp(1j*2*pi*nu_true*n*Tc);
                signal = alpha * phase;
                noise = sqrt(noise_var/2)*(randn(1,J)+1j*randn(1,J));
                received = signal + noise;
                
                Nfft = No * J;
                v = (0:Nfft-1) / (No*J*Tc);
                X = fft(received, Nfft);
                [~, idx_max] = max(abs(X));
                nu_est = v(idx_max);
                err_sum = err_sum + abs(nu_est - nu_true);
            end
            error(idx_snr, j_idx, no_idx) = err_sum / N_iter;
        end
    end
end

figure;
semilogy(SNR_dB, squeeze(error(:,1,1)), 'b-o'); hold on;
semilogy(SNR_dB, squeeze(error(:,2,1)), 'r-s');
semilogy(SNR_dB, squeeze(error(:,2,2)), 'g-d');
semilogy(SNR_dB, squeeze(error(:,2,3)), 'm-^');
xlabel('SNR (dB)'); ylabel('Avg Doppler Error (Hz)');
legend('J=5, No=10', 'J=10, No=10', 'J=10, No=100', 'J=10, No=1000');
grid on; title('Doppler Estimation Error');