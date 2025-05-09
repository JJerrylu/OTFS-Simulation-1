%% OTFS Simulation using DZT-OTFS with Multipath Channel, AWGN, and Pilot-Based Channel Estimation
clear all; close all; clc;

%% Parameters
L = 3;            % Number of delay bins
K = 4;            % Number of Doppler bins
N = L * K;        % Total number of time samples
SNR_dB = 20;      % SNR in dB for the AWGN channel

%% Generate delay-Doppler domain symbols (e.g., QPSK) for data transmission
% QPSK constellation: {1+j, 1-j, -1+j, -1-j}
constellation = [1+1j, 1-1j, -1+1j, -1-1j];
rng(123); % For reproducibility
Z_tx = constellation(randi([1, length(constellation)], K, L));

%% Transmit: IDZT to obtain time-domain signal X_tx from data grid
X_tx = zeros(N,1); % Preallocate time-domain signal

for n = 0:N-1
    l_idx = mod(n, L);         % Delay index: 0,...,L-1
    m = floor(n/L);            % Frame index (for Doppler phase)
    sum_val = 0;
    for k = 0:K-1
        % MATLAB indices start at 1, so add 1
        sum_val = sum_val + Z_tx(k+1, l_idx+1) * exp(1j*2*pi*m*(k/K));
    end
    X_tx(n+1) = (1/sqrt(K)) * sum_val;
end

%% Reshape transmitted time-domain signal to show delay-time structure
% Each column corresponds to a frame (m) and each row corresponds to delay index l.
X_tx_mat = reshape(X_tx, L, K);

%% Plot the transmitted delay-Doppler grid (data)
figure;
subplot(2,2,1);
bar3(real(Z_tx));
title('Tx Data: Delay-Doppler (Real)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Amplitude');

subplot(2,2,2);
bar3(imag(Z_tx));
title('Tx Data: Delay-Doppler (Imag)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Amplitude');

subplot(2,2,3);
stem(abs(X_tx));
title('Tx Data: Time-Domain Signal Magnitude (Serial)');
xlabel('n'); ylabel('|X[n]|');

subplot(2,2,4);
bar3(abs(X_tx_mat.'));
title('Tx Data: Delay-Time Domain Signal');
xlabel('Delay index (l)'); ylabel('Frame index (m)'); zlabel('Magnitude');

%% Define Multipath Channel Parameters
% We define a channel with 3 paths.
% Path 1: delay = 0, Doppler = 0, gain = 1
% Path 2: delay = 1, Doppler = +0.1, gain = 0.8
% Path 3: delay = 2, Doppler = -0.1, gain = 0.5
delays = [0, 1, 2];                   % in sample delays (corresponding to delay bins)
doppler_shifts = [0, 0.1, -0.1];        % normalized Doppler shifts (cycles per sample)
gains = [1, 0.8, 0.5];                % complex gains (assumed real here for simplicity)

% Also, for reference, define the true channel impulse response in delay-Doppler domain.
H_true = zeros(K, L);
H_true(1,1) = 1;        % Path 1 at (0,0)
H_true(2,2) = 0.8;      % Path 2 at (k=1, l=1) [mapping Doppler shift +0.1 to bin 2]
H_true(4,3) = 0.5;      % Path 3 at (k= -1, l=2) [for K=4, bin 4 represents a negative Doppler]
      
%% Simulate multipath propagation Data Signal in Time Domain
% Each path: circular shift by delay and apply Doppler phase modulation.
Y_channel = zeros(N,1);  % initialize channel output

n_vec = (0:N-1).';  % time index column vector
for p = 1:length(gains)
    % Circularly shift X_tx by delay (simulate delay effect)
    X_shifted = circshift(X_tx, delays(p));
    % Apply Doppler modulation: exp(j2*pi*doppler_shift*n)
    doppler_mod = exp(1j*2*pi*doppler_shifts(p)*n_vec);
    Y_channel = Y_channel + gains(p)*doppler_mod.*X_shifted;
end

%% Add AWGN to Data Channel Output
signal_power = mean(abs(Y_channel).^2);
SNR = 10^(SNR_dB/10);
noise_power = signal_power/SNR;
noise = sqrt(noise_power/2) * (randn(size(Y_channel)) + 1j*randn(size(Y_channel)));
Y_rx = Y_channel + noise;  % Received data signal (time domain)

%% Reshape Received Data Signal to Delay-Time Matrix
Y_rx_mat = reshape(Y_rx, L, K);

%% Receiver Processing for Data: Convert Y_rx back to Delay-Doppler domain via DZT
Z_rx = zeros(K,L);
for l = 0:L-1
    y_m = Y_rx_mat(l+1, :);  % 1 x K vector (one symbol per frame)
    Z_rx(:,l+1) = (1/sqrt(K)) * fft(y_m.');
end

%% Plot the Recovered Data Delay-Doppler Grid at Receiver
figure;
subplot(2,2,1);
bar3(real(Z_rx));
title('Rx Data: Delay-Doppler (Real)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Amplitude');

subplot(2,2,2);
bar3(imag(Z_rx));
title('Rx Data: Delay-Doppler (Imag)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Amplitude');

subplot(2,2,3);
stem(abs(Y_rx));
title('Rx Data: Time-Domain Signal Magnitude (Serial)');
xlabel('n'); ylabel('|Y[n]|');

subplot(2,2,4);
bar3(abs(Y_rx_mat'));
title('Rx Data: Delay-Time Domain Signal');
xlabel('Delay index (l)'); ylabel('Frame index (m)'); zlabel('Magnitude');

%% Pilot-Based Channel Estimation
% Create a pilot grid: all zeros except one pilot symbol at (k,l) = (0,0)
Z_pilot = zeros(K, L);
Z_pilot(1,1) = 1;   % Pilot placed at top-left (delay=0, Doppler=0)

% IDZT for the pilot grid to obtain time-domain pilot signal X_pilot
X_pilot = zeros(N,1);
for n = 0:N-1
    l_idx = mod(n, L);
    m = floor(n/L);
    sum_val = 0;
    for k = 0:K-1
        sum_val = sum_val + Z_pilot(k+1, l_idx+1)*exp(1j*2*pi*m*(k/K));
    end
    X_pilot(n+1) = (1/sqrt(K))*sum_val;
end

%% Transmit the Pilot Signal through the Multipath Channel
Y_pilot_channel = zeros(N,1);
for p = 1:length(gains)
    X_pilot_shifted = circshift(X_pilot, delays(p));
    doppler_mod = exp(1j*2*pi*doppler_shifts(p)*n_vec);
    Y_pilot_channel = Y_pilot_channel + gains(p)*doppler_mod.*X_pilot_shifted;
end

% Add AWGN to the pilot channel output
signal_power_pilot = mean(abs(Y_pilot_channel).^2);
noise_power_pilot = signal_power_pilot/SNR;
noise_pilot = sqrt(noise_power_pilot/2) * (randn(size(Y_pilot_channel)) + 1j*randn(size(Y_pilot_channel)));
Y_pilot = Y_pilot_channel + noise_pilot;

%% Receiver Processing for Pilot: Convert received pilot time-domain signal to Delay-Doppler via DZT
Y_pilot_mat = reshape(Y_pilot, L, K);
Z_pilot_rx = zeros(K,L);
for l = 0:L-1
    y_m = Y_pilot_mat(l+1, :);
    Z_pilot_rx(:,l+1) = (1/sqrt(K)) * fft(y_m.');
end

%% The estimated channel impulse response in delay-Doppler domain is given by the recovered pilot grid.
H_est = Z_pilot_rx;  % (For a properly designed pilot, H_est corresponds to the channel response (possibly shifted))

%% Plot the Estimated Channel Impulse Response (3D Bar Chart)
figure;
subplot(1,2,1);
bar3(abs(H_true));
title('Real Channel Impulse Response (Magnitude)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('|H|');

subplot(1,2,2);
bar3(abs(H_est));
title('Estimated Channel Impulse Response (Magnitude) - Image');
xlabel('Delay index (l)'); ylabel('Doppler index (k)');

%% Coarse channel equalization

Z_X_est = Z_rx./H_est;

%% Plot the equalized received data (3D Bar Chart)
figure;
subplot(2,2,1);
bar3(real(Z_tx));
title('Tx Data: Delay-Doppler (Real)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Amplitude');

subplot(2,2,2);
bar3(imag(Z_tx));
title('Tx Data: Delay-Doppler (Imag)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)'); zlabel('Amplitude');


subplot(2,2,3);
bar3(real(Z_X_est));
title('Estimated data: (Delay-Doppler (Real)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)');

subplot(2,2,4);
bar3(imag(Z_X_est));
title('Estimated data: (Delay-Doppler (image)');
xlabel('Delay index (l)'); ylabel('Doppler index (k)');


%% Display Key Matrices for Reference
disp('Transmitted Data Delay-Doppler Grid (Z_tx):');
disp(Z_tx);
disp('Received Data Delay-Doppler Grid (Z_rx):');
disp(Z_rx);
disp('True Channel Impulse Response in Delay-Doppler Domain (H_true):');
disp(H_true);
disp('Estimated Channel Impulse Response in Delay-Doppler Domain (H_est):');
disp(H_est);
