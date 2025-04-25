% Parameters
numSymbols = 100;         % Number of data symbols
sps = 8;                  % Samples per symbol (oversampling factor)
span = 10;                % Filter span in symbols
M = 4;                    % Modulation order (QPSK)
data = randi([0 M-1], numSymbols, 1);
modData = pskmod(data, M, pi/M);  % QPSK modulation
upsampled = zeros(numSymbols * sps, 1);
upsampled(1:sps:end) = modData;  % Insert symbols with zeros between
% Time vector for filter
t = (-span*sps/2 : span*sps/2) / sps;  % In symbol durations

% Sinc pulse (normalized)
sincPulse = sinc(t);  % sinc(x) = sin(pi*x)/(pi*x)

% Optional: apply a window to control sidelobes (e.g., Hamming)
window = hamming(length(sincPulse))';
sincPulse = sincPulse .* window;
txSignal = conv(upsampled, sincPulse, 'same');
% Time-domain plot
figure;
plot(real(txSignal));
title('Time-domain Transmit Signal (Sinc Pulse Shaped)');
xlabel('Sample Index'); ylabel('Amplitude');

% Frequency-domain plot
figure;
pwelch(txSignal, [], [], [], sps*1e3, 'centered');
title('Frequency-domain Spectrum');
