%% OTFS Simulation using DZT-OTFS with AWGN Channel
clear all; close all; clc;

%% Parameters
L = 3;            % Number of delay bins
K = 4;            % Number of Doppler bins
N = L * K;        % Total number of time samples
SNR_dB = 10;      % SNR in dB for the AWGN channel

%% Generate delay-Doppler domain symbols (e.g., QPSK)
% Here we use QPSK symbols from the set {1+j, 1-j, -1+j, -1-j}
constellation = [1+1j, 1-1j, -1+1j, -1-1j];
rng(123); % for reproducibility
Z_tx = constellation(randi([1, length(constellation)], K, L));

%% Transmit: IDZT to obtain time-domain signal X[n]
X_tx = zeros(N,1); % preallocate time-domain signal

for n = 0:N-1
    l_idx = mod(n, L);         % delay index (0,...,L-1)
    m = floor(n/L);            % frame (or time index) corresponding to Doppler phase
    sum_val = 0;
    for k = 0:K-1
        % Note: MATLAB indices start at 1, so add 1
        sum_val = sum_val + Z_tx(k+1, l_idx+1) * exp(1j*2*pi*m*(k/K));
    end
    X_tx(n+1) = (1/sqrt(K)) * sum_val;
end

%% Reshape transmitted time-domain signal to show delay-time structure
% Each column corresponds to a frame (m) and each row corresponds to delay index l.
X_tx_mat = reshape(X_tx, L, K);

%% Plot the transmitted delay-Doppler grid
figure;
subplot(2,2,1);
bar3(real(Z_tx));
title('Tx Delay-Doppler (Real part)');
xlabel('Delay index (l)');
ylabel('Doppler index (k)');

subplot(2,2,2);
bar3(imag(Z_tx));
title('Tx Delay-Doppler (Imag part)');
xlabel('Delay index (l)');
ylabel('Doppler index (k)');

subplot(2,2,3);
stem(abs(X_tx_mat(:)));
title('Tx Time-Domain Signal Magnitude (Serial)');
xlabel('n'); ylabel('|X[n]|');

subplot(2,2,4);
bar3(abs(X_tx_mat.'));
title('Tx Delay-Time Domain Signal');
xlabel('Delay index l'); ylabel('Frame index m');

%% Pass the transmitted signal through an AWGN channel
% Compute noise power from SNR_dB
signal_power = mean(abs(X_tx).^2);
SNR = 10^(SNR_dB/10);
noise_power = signal_power/SNR;
noise = sqrt(noise_power/2) * (randn(size(X_tx)) + 1j*randn(size(X_tx)));
Y_rx = X_tx + noise;  % Received time-domain signal

%% Reshape received time-domain signal to delay-time matrix
Y_rx_mat = reshape(Y_rx, L, K);

%% Receiver Processing: Convert Y_rx back to Delay-Doppler domain via DZT
% For each delay bin l and Doppler index k, perform a K-point DFT along the frame axis.
Z_rx = zeros(K,L);
for l = 0:L-1
    % Extract the m-indexed samples for this delay index l:
    y_m = Y_rx_mat(l+1, :);  % This is a 1 x K vector (one symbol per frame)
    % Perform a K-point DFT with proper scaling:
    Z_rx(:,l+1) = (1/sqrt(K)) * fft(y_m.');   % Transpose to get a column vector
end

%% Plot the recovered delay-Doppler grid at receiver
figure;
subplot(2,2,1);
bar3(real(Z_rx));
title('Rx Delay-Doppler (Real part)');
xlabel('Delay index (l)');
ylabel('Doppler index (k)');

subplot(2,2,2);
bar3(imag(Z_rx));
title('Rx Delay-Doppler (Imag part)');
xlabel('Delay index (l)');
ylabel('Doppler index (k)');

%% Also plot the received time-domain signal (delay-time domain)
subplot(2,2,3);
stem(abs(Y_rx)); 
title('Rx Time-Domain Signal Magnitude (Serial)');
xlabel('n'); ylabel('|Y[n]|');

subplot(2,2,4);
bar3(abs(Y_rx_mat'));
title('Rx Delay-Time Domain Signal');
xlabel('Delay index l'); ylabel('Frame index m');

%% Display constellation error (if desired)
disp('Transmitted Delay-Doppler Grid (Z_tx):');
disp(Z_tx);
disp('Received Delay-Doppler Grid (Z_rx):');
disp(Z_rx);
