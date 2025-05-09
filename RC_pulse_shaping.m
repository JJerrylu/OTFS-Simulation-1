% ------------------------------
% Pulse shaping: sinc-RC pulse
% ------------------------------
function plot_sinc_rc_pulse()
    Ts = 1/60e3;          % Symbol duration
    osr = 10;             % Oversampling rate
    t_range = -4*Ts : Ts/osr : 4*Ts;  % Time vector (oversampled)

    gammas = [0, 0.25, 0.5];
    colors = ['b', 'r', 'g'];
    figure; hold on; grid on;
    for i = 1:length(gammas)
        gamma = gammas(i);
        p = pulse_shape(t_range, Ts, gamma);
        plot(t_range * 1e6, p, colors(i), 'LineWidth', 2);  % 時間以 us 為單位
    end
    xlabel('Time (    xlabel('Time (\x03bcs)'); ylabel('Amplitude');
    title('Sinc-RC Pulse with Different Roll-Off Factors');
    legend('\gamma = 0', '\gamma = 0.25', '\gamma = 0.5');
end

% ------------------------------
% Pulse shape function
% ------------------------------
function p = pulse_shape(t, Ts, gamma)
    % Avoid singularity at denominator
    denom = 1 - (2 * gamma * t / Ts).^2;
    denom(abs(denom) < 1e-8) = 1e-8;

    sinc_term = sinc(t / Ts);  % MATLAB sinc = sin(pi*x)/(pi*x)
    cosine_term = cos(pi * gamma * t / Ts);

    p = sinc_term .* cosine_term ./ denom;
    p = p / norm(p);  % Normalize energy if needed
end
