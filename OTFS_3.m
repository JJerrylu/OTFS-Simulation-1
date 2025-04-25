%% OTFS Simulation with Embedded Pilot, Iterative Channel Estimation, and Equalization
clear all; close all; clc;

%% Parameters
L = 3;            % Number of delay bins
K = 4;            % Number of Doppler bins
N = L * K;        % Total number of time samples
SNR_dB = 10;      % SNR in dB for the AWGN channel

%% Generate Data Grid and Embed Pilot
% QPSK constellation: {1+j, 1-j, -1+j, -1-j}
constellation = [1+1j, 1-1j, -1+1j, -1-1j];
rng(123); % For reproducibility
% Generate random data for all cells
Z_data = constellation(randi([1, length(constellation)], K, L));
% Embed pilot: For example, reserve cell (1,1) (i.e. delay=0, doppler=0) for pilot
pilot_val = 1;
Z_data(1,1) = pilot_val;
% This is the transmitted delay-Doppler grid (data + embedded pilot)
Z_tx = Z_data;

%% Transmit: Map Delay-Doppler Grid to Time Domain via IDZT
X_tx = zeros(N,1); % Preallocate time-domain signal
for n = 0:N-1
    l_idx = mod(n, L);         % delay index (0,...,L-1)
    m = floor(n/L);            % frame index (for Doppler phase)
    sum_val = 0;
    for k = 0:K-1
        % Note: MATLAB indices start at 1
        sum_val = sum_val + Z_tx(k+1, l_idx+1) * exp(1j*2*pi*m*(k/K));
    end
    X_tx(n+1) = (1/sqrt(K))*sum_val;
end
X_tx_mat = reshape(X_tx, L, K);

%% Plot Transmitted Signals (Data and Pilot Embedded)
figure;
subplot(2,2,1);
bar3(real(Z_tx));
title('Tx Delay-Doppler Grid (Real)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Amplitude');

subplot(2,2,2);
bar3(imag(Z_tx));
title('Tx Delay-Doppler Grid (Imag)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Amplitude');

subplot(2,2,3);
stem(abs(X_tx));
title('Tx Time-Domain Signal (Serial)');
xlabel('n'); ylabel('|X[n]|');

subplot(2,2,4);
bar3(abs(X_tx_mat.'));
title('Tx Delay-Time Domain Signal');
xlabel('Delay index (l)'); ylabel('Frame index (m)'); zlabel('Magnitude');

%% Define Multipath Channel Parameters
% Assume a channel with three paths:
%  Path 1: delay = 0 samples, Doppler = 0, gain = 1
%  Path 2: delay = 1 sample, Doppler = +0.1, gain = 0.8
%  Path 3: delay = 2 samples, Doppler = -0.1, gain = 0.5
delays = [0, 1, 2];                   % in sample delays (corresponding to delay bins)
doppler_shifts = [0, 0.1, -0.1];        % normalized Doppler shifts (cycles per sample)
gains = [1, 0.8, 0.5];                % complex gains (here assumed real)

% For reference, we define the true channel impulse response (in delay-Doppler grid)
H_true = zeros(K, L);
H_true(1,1) = 1;        % Path 1 at (doppler=0, delay=0)
H_true(2,2) = 0.8;      % Path 2: here we map Doppler +0.1 to bin 2 and delay=1 to l=1
H_true(4,3) = 0.5;      % Path 3: for K=4, bin 4 represents a negative doppler (e.g. -1) and delay=2 to l=2

%% Channel Simulation in Time Domain (Pilot and Data Transmitted Together)
Y_channel = zeros(N,1);  % initialize channel output
n_vec = (0:N-1).';       % time indices as a column vector

for p = 1:length(gains)
    % Apply circular delay shift and Doppler modulation for each path:
    X_shifted = circshift(X_tx, delays(p));
    doppler_mod = exp(1j*2*pi*doppler_shifts(p)*n_vec);
    Y_channel = Y_channel + gains(p) * doppler_mod .* X_shifted;
end

% Add AWGN noise
signal_power = mean(abs(Y_channel).^2);
SNR = 10^(SNR_dB/10);
noise_power = signal_power/SNR;
noise = sqrt(noise_power/2)*(randn(size(Y_channel)) + 1j*randn(size(Y_channel)));
Y_rx = Y_channel + noise;  % Received time-domain signal

%% Reshape Received Time-Domain Signal to Delay-Time Matrix
Y_rx_mat = reshape(Y_rx, L, K);

%% Receiver Processing: Map Received Time-Domain Signal to Delay-Doppler Domain via DZT
Z_rx_total = zeros(K, L);
for l = 0:L-1
    y_m = Y_rx_mat(l+1, :);  % one row corresponding to a given delay bin over frames
    Z_rx_total(:,l+1) = (1/sqrt(K)) * fft(y_m.');  % perform K-point DFT along frame dimension
end

%% Embedded Pilot Based Channel Estimation
% In an embedded pilot scheme, the known pilot (here at cell (1,1)) is used to
% extract a rough channel impulse response. In practical systems a larger pilot block
% is used so that the delay-Doppler spread can be estimated.
% Here we assume that the channel’s effect on the pilot is spread over the grid.
% A simple (and rough) estimator is to take the entire received grid (Z_rx_total)
% and divide by the known pilot value (assuming interference is low).
H_est = Z_rx_total / pilot_val;  

%% Plot the Estimated Channel Impulse Response
figure;
subplot(1,2,1);
bar3(abs(H_est));
title('Estimated Channel Impulse Response (Magnitude)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('|H|');
subplot(1,2,2);
bar3(abs(H_est));
title('Estimated Channel Impulse Response (Image)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)');

%% Iterative Detection / Peeling of Channel Paths
% Here we perform a simple iterative procedure that, for a known number of paths,
% finds the strongest tap in H_est, records its parameters, and “cancels” it from
% the estimate before moving on.
H_iter = H_est;   % working copy of the estimated channel response
num_paths = 3;    % assume we know there are three dominant paths
estimated_paths = [];  % to store [gain, doppler_bin, delay_bin]

for iter = 1:num_paths
    [max_val, idx] = max(H_iter(:));
    [k_idx, l_idx] = ind2sub(size(H_iter), idx);  % indices in H_iter (1-indexed)
    % Convert indices to (doppler, delay) bin numbers (0-indexed)
    est_gain = max_val;
    est_doppler_bin = k_idx - 1;
    est_delay_bin = l_idx - 1;
    estimated_paths = [estimated_paths; est_gain, est_doppler_bin, est_delay_bin];
    % Cancel this path (set the corresponding cell to zero)
    H_iter(k_idx, l_idx) = 0;
end
disp('Estimated Channel Paths [gain, doppler bin, delay bin]:');
disp(estimated_paths);

%% Data Equalization Using the Estimated Channel
% A simple equalizer (for illustration) divides the received delay-Doppler signal by
% the estimated channel response at each cell (skipping the pilot cell).
Z_eq = Z_rx_total;
for k = 1:K
    for l = 1:L
        if (k==1 && l==1)
            continue;  % leave the pilot cell unchanged
        end
        if abs(H_est(k,l)) > 1e-3
            Z_eq(k,l) = Z_rx_total(k,l) / H_est(k,l);
        end
    end
end

%% Decision: Detect Data Symbols by Nearest Constellation Point
Z_detect = Z_eq;
for k = 1:K
    for l = 1:L
        if (k==1 && l==1)
            continue;  % pilot is known
        end
        distances = abs(Z_eq(k,l) - constellation);
        [~, min_idx] = min(distances);
        Z_detect(k,l) = constellation(min_idx);
    end
end

%% Display Results
disp('Transmitted Delay-Doppler Grid (with embedded pilot):');
disp(Z_tx);
disp('Equalized and Detected Delay-Doppler Grid:');
disp(Z_detect);

%% Plot the Detected Data Delay-Doppler Grid
figure;
subplot(1,2,1);
bar3(abs(Z_detect));
title('Detected Data Delay-Doppler Grid (Magnitude)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Magnitude');
subplot(1,2,2);
bar3(abs(Z_detect)); 
title('Detected Data Delay-Doppler Grid (Image)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)');
