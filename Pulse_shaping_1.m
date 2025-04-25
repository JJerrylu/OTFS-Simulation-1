% Parameters
sps = 8;                 % Samples per symbol
span = 10;               % Symbols
N = span * sps + 1;      % Total number of samples
t = (-span/2 : 1/sps : span/2);  % Time vector

% Ideal sinc pulse
sincPulse = sinc(t);

% Define windows
rectWindow = ones(size(sincPulse));              % Rectangular (no taper)
hammingWindow = hamming(length(sincPulse))';     % Hamming window
hannWindow = hann(length(sincPulse))';           % Hann window

% Apply windows
sincRect = sincPulse .* rectWindow;
sincHamming = sincPulse .* hammingWindow;
sincHann = sincPulse .* hannWindow;

% Plot Time Domain
figure;
subplot(3,1,1); plot(t, sincRect); title('Sinc with Rectangular Window'); xlabel('Time'); ylabel('Amplitude');
subplot(3,1,2); plot(t, sincHamming); title('Sinc with Hamming Window'); xlabel('Time'); ylabel('Amplitude');
subplot(3,1,3); plot(t, sincHann); title('Sinc with Hann Window'); xlabel('Time'); ylabel('Amplitude');

% Plot Frequency Domain
figure;
[H1, f] = freqz(sincRect, 1, 1024, sps*1e3);
[H2, ~] = freqz(sincHamming, 1, 1024, sps*1e3);
[H3, ~] = freqz(sincHann, 1, 1024, sps*1e3);

plot(f, 20*log10(abs(H1)), 'r', 'DisplayName','Rectangular');
hold on;
plot(f, 20*log10(abs(H2)), 'b', 'DisplayName','Hamming');
plot(f, 20*log10(abs(H3)), 'g', 'DisplayName','Hann');
legend;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response of Windowed Sinc Filters');
grid on;
