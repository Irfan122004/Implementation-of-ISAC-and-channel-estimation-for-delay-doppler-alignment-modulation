clear; clc; close all;

B = 100e6; Ts = 1/B;
N_vals = [5000, 10000, 100000];
d_tau = -50:50;
d_nu = -5000:100:5000;

figure;
for n_idx = 1:length(N_vals)
    N = N_vals(n_idx);
    % Doppler-cut AF
    chi_doppler = zeros(size(d_nu));
    for idx = 1:length(d_nu)
        if d_nu(idx) == 0
            chi_doppler(idx) = 1;
        else
            chi_doppler(idx) = abs(sin(pi*d_nu(idx)*N*Ts) / (N*sin(pi*d_nu(idx)*Ts)));
        end
    end
    subplot(3,2,2*n_idx-1);
    plot(d_nu, 20*log10(chi_doppler), 'b-', 'LineWidth', 1.5);
    xlabel('Doppler difference (Hz)'); ylabel('Normalized AF (dB)');
    title(sprintf('Doppler-cut AF, N=%d', N));
    grid on; ylim([-40 5]);
    
    % Delay-cut AF (single-path case -> zero sidelobes)
    chi_delay = zeros(size(d_tau));
    chi_delay(d_tau==0) = 1;
    subplot(3,2,2*n_idx);
    stem(d_tau, chi_delay, 'filled');
    xlabel('Delay difference (samples)'); ylabel('Normalized AF');
    title(sprintf('Delay-cut AF, N=%d', N));
    grid on; ylim([0 1.2]);
end