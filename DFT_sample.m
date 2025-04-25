% Parameters
Fs = 1000;        % Sampling frequency (Hz)
Ts = 1/Fs;        % Sampling period (s)
N = 128;            % Number of samples (DFT size)
t = (0:N-1)*Ts;   % Time vector (duration = N*Ts)

% Time-domain signal (e.g., simple sinusoid)
f_signal = 2500;                     % Signal frequency (Hz)
x = cos(2*pi*f_signal*t);          % Real-valued time-domain signal

% Compute DFT
X = fft(x);                         % N-point DFT
f = (0:N-1)*(Fs/N);                 % Frequency vector (raw indices)

% Shift zero frequency to center (for plotting with negative freqs)
X_shifted = fftshift(X);
f_shifted = (-N/2:N/2-1)*(Fs/N);   % Frequency vector centered at 0 Hz

% Display time and frequency resolution
T_total = N * Ts;                 % Total duration
delta_f = 1 / T_total;           % Frequency resolution

fprintf('Total duration = %.4f s\n', T_total);
fprintf('Frequency resolution = %.2f Hz\n', delta_f);

% Plot
subplot(2,1,1);
plot(t, x);
title('Time Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(f_shifted, abs(X_shifted));
title('Frequency Domain (Magnitude Spectrum)');
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
grid on;
