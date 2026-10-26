clear; clc; close all;

M = 64; P = 100; L = 5; J_vals = [1, 5, 10];
SNR_dB = -10:5:20;
N_iter = 100;

NMSE_OMP = zeros(length(SNR_dB), 1);
NMSE_SOMP = zeros(length(SNR_dB), length(J_vals));
NMSE_ASOMP = zeros(length(SNR_dB), length(J_vals));

for idx_snr = 1:length(SNR_dB)
    SNR = 10^(SNR_dB(idx_snr)/10);
    noise_var = 1/SNR;
    
    for iter = 1:N_iter
        H_true = zeros(M, P);
        path_delays = randi(P, L, 1);
        for l = 1:L
            theta = (rand-0.5)*pi;
            a = exp(1j*pi*sin(theta)*(0:M-1)')./sqrt(M);
            H_true(:, path_delays(l)) = a * (randn + 1j*randn);
        end
        h_true = H_true(:);
        
        Np = 100;
        Phi = randn(Np, M*P) + 1j*randn(Np, M*P);
        Phi = Phi ./ sqrt(sum(abs(Phi).^2,1));
        
        for j_idx = 1:length(J_vals)
            J = J_vals(j_idx);
            Y = zeros(Np, J);
            for k = 1:J
                y = Phi * h_true + sqrt(noise_var/2)*(randn(Np,1)+1j*randn(Np,1));
                Y(:,k) = y;
            end
            
            h_est_omp = omp_est(Phi, Y(:,1), L);
            NMSE_OMP(idx_snr) = NMSE_OMP(idx_snr) + norm(h_est_omp - h_true)^2 / norm(h_true)^2;
            
            h_est_somp = somp(Phi, Y, L);
            NMSE_SOMP(idx_snr, j_idx) = NMSE_SOMP(idx_snr, j_idx) + norm(h_est_somp - h_true)^2 / norm(h_true)^2;
            
            h_est_asomp = asomp_sr(Phi, Y, L, 8, 8);
            NMSE_ASOMP(idx_snr, j_idx) = NMSE_ASOMP(idx_snr, j_idx) + norm(h_est_asomp - h_true)^2 / norm(h_true)^2;
        end
    end
end

NMSE_OMP = NMSE_OMP / N_iter;
NMSE_SOMP = NMSE_SOMP / N_iter;
NMSE_ASOMP = NMSE_ASOMP / N_iter;

figure;
semilogy(SNR_dB, NMSE_OMP, 'b-o', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB, NMSE_SOMP(:,1), 'r--s', 'LineWidth', 1.5);
semilogy(SNR_dB, NMSE_SOMP(:,2), 'g--^', 'LineWidth', 1.5);
semilogy(SNR_dB, NMSE_ASOMP(:,2), 'm-d', 'LineWidth', 1.5);
xlabel('SNR (dB)'); ylabel('NMSE');
legend('OMP', 'SOMP J=1', 'SOMP J=5', 'ASOMP-SR J=5');
grid on; title('NMSE Performance (On-Grid)');


function h_est = omp_est(Phi, y, L)
    [N, MP] = size(Phi);
    residual = y;
    support = [];
    for iter = 1:L
        corr = abs(Phi' * residual);
        [~, idx] = max(corr);
        support = [support, idx];
        Phi_s = Phi(:, support);
        h_temp = Phi_s \ y;
        residual = y - Phi_s * h_temp;
    end
    h_est = zeros(MP,1);
    h_est(support) = h_temp;
end

function h_est = somp(Phi, Y, L)
    [N, MP] = size(Phi);
    [~, J] = size(Y);
    residual = Y;
    support = [];
    for iter = 1:L
        corr_sum = zeros(MP,1);
        for k = 1:J
            corr_sum = corr_sum + abs(Phi' * residual(:,k)).^2;
        end
        [~, idx] = max(corr_sum);
        support = [support, idx];
        Phi_s = Phi(:, support);
        h_temp = Phi_s \ Y;
        residual = Y - Phi_s * h_temp;
    end
    h_est = zeros(MP, J);
    h_est(support, :) = h_temp;
    h_est = h_est(:,1);
end

function h_est = asomp_sr(Phi, Y, L, Vr, Vp)
    [N, MP] = size(Phi);
    [~, J] = size(Y);
    residual = Y;
    support = [];
    for iter = 1:L
        corr_sum = zeros(MP,1);
        for k = 1:J
            corr_sum = corr_sum + abs(Phi' * residual(:,k)).^2;
        end
        [~, idx] = max(corr_sum);
        support = [support, idx];
        Phi_s = Phi(:, support);
        h_temp = Phi_s \ Y;
        residual = Y - Phi_s * h_temp;
    end
    refined_support = [];
    for idx = support
        p = floor((idx-1)/MP);  % simplified delay index
        r = mod(idx-1, MP);      % simplified angular index
        for dr = -Vr/2:Vr/2-1
            for dp = -Vp/2:Vp/2-1
                new_idx = (p+dp)*MP + (r+dr) + 1;
                if new_idx >= 1 && new_idx <= MP
                    refined_support = [refined_support, new_idx];
                end
            end
        end
    end
    refined_support = unique(refined_support);
    Phi_s = Phi(:, refined_support);
    h_temp = Phi_s \ Y;
    h_est = zeros(MP, J);
    h_est(refined_support, :) = h_temp;
    h_est = h_est(:,1);
end