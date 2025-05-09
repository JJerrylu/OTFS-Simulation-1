% --- Parameters ---
Ts = 1;            % Symbol duration
gamma = 0;       % Roll-off factor
span = 6;          % Pulse span in symbol durations
L = 10;             % Oversampling factor (samples per Ts)
num_symbols = 100;  % Number of symbols

% --- Generate random BPSK symbols ---
%symbols = 2 * randi([0 1], 1, num_symbols) - 1;   % BPSK: ±1
symbols =ones(1,num_symbols);
% --- Generate raised cosine pulse ---
g = raised_cosine_pulse(Ts, gamma, span, L);

% --- Upsample symbol sequence ---
upsampled = zeros(1, num_symbols * L);
upsampled(1:L:end) = symbols;  % Insert zeros between symbols

% --- Perform pulse shaping by convolution ---
tx_signal = conv(upsampled, g, 'full');  % shaped baseband signal

% --- Time axis for plotting ---
T_total = length(tx_signal);
t_axis = (0:T_total-1) * (Ts / L);  % time in seconds

% --- Plot results ---
figure;
plot(t_axis, tx_signal, 'b'); grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Pulse-Shaped Baseband Signal');

function g = raised_cosine_pulse(Ts, gamma, span, L)
% Ts     - Symbol duration
% gamma  - Roll-off factor (0 <= gamma <= 1)
% span   - Pulse span in symbol durations (e.g., 6 means [-3T, 3T])
% L      - Oversampling factor (samples per Ts)

% Total pulse duration
T_total = span * Ts;                 
t = -T_total/2 : Ts/L : T_total/2;   % time axis

% Initialize g(t)
g = zeros(size(t));

for i = 1:length(t)
    ti = t(i);
    
    % Handle t == 0 (singularity)
    if abs(ti) < 1e-8
        g(i) = 1;
    % Handle denominator zero (sinc singularity)
    elseif abs(abs(2 * gamma * ti / Ts) - 1) < 1e-8
        g(i) = (pi/4) * sinc(1/(2 * gamma));
    else
        numerator = sin(pi * ti / Ts) .* cos(gamma * pi * ti / Ts);
        denominator = (pi * ti / Ts) .* (1 - (2 * gamma * ti / Ts)^2);
        g(i) = numerator / denominator;
    end
end

% Normalize energy (optional)
g = g / norm(g);

% Plot (optional)
% figure; plot(t, g);
% xlabel('Time'); ylabel('g(t)'); title('Raised Cosine Pulse');
end

