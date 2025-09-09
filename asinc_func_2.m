clc; clear all; close all;

N = 32;                  % Doppler grid size
L = 100;                 % oversampling factor
N_os = (N-1)*L+1;        % oversampled points

frac = [3.5, -4.7];   % fractional Doppler shifts

% oversampled Doppler axis
k_os = linspace(-N/2, N/2-1, N_os);

figure; hold on;
legend_entries = {};     % collect legend text

for ii = 1:length(frac)
    kernel_os = zeros(1, N_os);
    for idx = 1:N_os
        k = k_os(idx);
        % Dirichlet kernel (asinc)
        num = exp(-1i*2*pi*(k - frac(ii))) - 1;
        den = exp(-1i*2*pi*(k - frac(ii))/N) - 1;
        if abs(den) < 1e-12
            kernel_os(idx) = 1;   % limit case (normalized)
        else
            kernel_os(idx) = num ./ den / N;
        end
    end

    % plot oversampled curve
    plot(k_os, abs(kernel_os), 'LineWidth', 1.5);

    % plot discrete samples
    stem(-N/2:N/2-1, abs(kernel_os(1:L:end)), 'filled');

    % save legend entry
    legend_entries{end+1} = sprintf("frac = %.2f", frac(ii));
end

grid on;
title("Fractional Doppler Asinc Kernel (Oversampled)");
xlabel("k (oversampled Doppler index)");
ylabel("|Kernel|");
legend(legend_entries, 'Location', 'best');
