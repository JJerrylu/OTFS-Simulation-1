clc; clear all; close all;

N = 32;                  % Doppler grid size
L = 100;                 % oversampling factor
N_os = (N-1)*L+1;        % oversampled points

frac = [8.5, -7.7, 0];   % fractional Doppler shifts
color1 = ['k', 'b','c'];
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
    plot(k_os, abs(kernel_os), 'LineWidth', 2,'Color',color1(ii));
    % save legend entry
    legend_entries{end+1} = sprintf("frac = %.1f", frac(ii));

    % plot discrete samples
    stem(-N/2:N/2-1, abs(kernel_os(1:L:end)),'r','LineWidth',1.5);
    legend_entries{end+1} = sprintf("Integer Sampled Point");
end

grid on;
title = title("Fractional Doppler Asinc Kernel (Oversampled)");
title.FontSize = 20;
title.FontWeight = "bold";
title.Interpreter = "latex";
xlabel = xlabel("k (oversampled Doppler index)");
xlabel.FontSize = 20;
xlabel.Interpreter = "latex";
ylabel = ylabel("|h|");
ylabel.FontSize = 20;
ylabel.Interpreter = "latex";
lgd = legend(legend_entries);
lgd.FontSize = 20;
lgd.FontWeight = "bold";
lgd.Box = "on";
lgd.Location = "best";
lgd.Interpreter = "latex";
