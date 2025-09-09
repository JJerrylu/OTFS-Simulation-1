clc, clear all, close all
M = 512;
CP_length = 64;
N = 128;
MN = M*N;
frame_length = MN+CP_length*N;
delta_f = 15e3;            % 15kHz
T = 1/delta_f;             % Block Duration
fc = 4e9;                  % 4GHz
delay_tap = [0, 1, 2, 3, 4, 5, 8, 13, 19];
%delay_tap = [0.1, 2.2, 4.3, 6.4, 8.5, 10.6, 12.7, 14.8, 16.9];
%int_delay_inx = round(delay_tap);
%frac_delay_inx = delay_tap - int_delay_inx;
tap = length(delay_tap);

%% Delay and channel coefficient
delay = delay_tap/(M*delta_f);
delay_pow_prof = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];
lin_pow_prof = 10.^(delay_pow_prof/10);
lin_pow_prof = lin_pow_prof/sum(lin_pow_prof);
chan_coef = sqrt(lin_pow_prof).*(sqrt(1/2) * (randn(1,tap)+1i*randn(1,tap)));

%% Doppler shift
max_speed = 500;           % 500 km/h
v = max_speed*(1000/3600);
max_doppler = (v*fc)/(299792458);
doppler = max_doppler*cos(2*pi*rand(1,tap));
doppler_inx = doppler*(N*T); % Doppler taps using Jake's spectrum
%doppler_inx = [1, 2, -4, -3, 5, 6, -10, 13, -14];
%doppler = doppler_inx/(N*T);
int_doppler_inx = round(doppler_inx);                  % Integer Doppler
frac_doppler_inx = doppler_inx-int_doppler_inx;        % Fractional Doppler

X_DD     = zeros(M, N);
pilot_power_dB = [10 20 30 40];       %dB
noise = sqrt(1/2)*(randn(MN, 1)+1i*randn(MN,1));    % noise power = 0dB

pilot_power    = 10^(pilot_power_dB/10);
X_DD(M/2, N/2)  = sqrt(pilot_power);
X_DD_vec = reshape(X_DD, MN, 1);
%s        = ifft(X_DD, N, 2); % delay time domain transmitted symbol
F_N     = dftmtx(N);
F_N     = F_N/norm(F_N);
F_N_H   = F_N';
s = X_DD*F_N_H;
s_vec    = reshape(s, MN, 1);     % time domain transmitted symbol
H   = zeros(MN, MN);

doppler_mtx = zeros(M,M);

for n = 1:N
    H_n = zeros(M, M);
    for p = 1:tap 
        delay_mtx = eye(M);
        for i = 1:M
            doppler_mtx(i,i) = exp(-1i*2*pi*doppler_inx(p)/frame_length*((M+CP_length)*(n-1)- delay_tap(p)+i-1));
            delay_mtx(:, i) = circshift(delay_mtx(:, i), delay_tap(p));
        end
        H_n = H_n + chan_coef(p)*doppler_mtx*delay_mtx;
    end
    
        H((n-1)*M+1:n*M, (n-1)*M+1:n*M) = H_n;
end
r_vec = H * s_vec + noise ;              % Vectorized Time domain IO relation
r     = reshape(r_vec, M, N);    % Delay-Time domain recieved symbol
r_vec_noise_free = H * s_vec;              % Vectorized Time domain IO relation
r_noise_free     = reshape(r_vec_noise_free, M, N);    % Delay-Time domain recieved symbol
F_NM    = kron(F_N, eye(M));
F_NM_H  = F_NM';

H_eq     = F_NM * H * F_NM_H;  % DD domain channel matrix
Y_DD_vec = H_eq * X_DD_vec + F_NM*noise;   % Vectorized DD domain IO relation
Y_DD = reshape(Y_DD_vec, M , N); % Received DD grid
%Y_DD_verify = fft(r, N, 2);
Y_DD_verify = r*F_N;
Y_DD_noise_free = r_noise_free*F_N;
%% Plot
%figure; bar3(abs(X_DD)); title('DD domain transmitted pilot X\_DD (abs)'); xlabel('doppler'); ylabel('delay');
%figure; bar3(abs(H_n));  title('Delay domain channel matrix of n\_th OFDM symbol  H\_n (abs)'); xlabel('doppler'); ylabel('delay');

figure;
imagesc(0:M-1, 0:N-1, abs(X_DD).^2);
set(gca,'XDir','normal'); 
title('DD domain transmitted symbol Heat Map');
colormap(jet); colorbar;

figure;
imagesc(0:M-1, 0:N-1, abs(s).^2);
set(gca,'XDir','normal'); 
title('Delay_time domain transmitted symbol Heat Map');
colormap(jet); colorbar;

figure;
imagesc(0:MN-1, 0:MN-1, abs(H_n).^2);
set(gca,'XDir','normal'); 
title('Delay domain channel matrix of n\_th OFDM symbol Heat Map');
colormap(jet); colorbar;

figure;
imagesc(0:MN-1, 0:MN-1, abs(H).^2);
set(gca,'XDir','normal'); 
title('Time domain Channel Magnitude Heat Map');
colormap(jet); colorbar;

figure;
imagesc(0:M-1, 0:N-1, abs(r).^2);
set(gca,'XDir','normal'); 
title('Delay-Time domain Received Symbol Heat Map');
colormap(jet); colorbar;

figure;
imagesc(0:M-1, 0:N-1, abs(r_noise_free).^2);
set(gca,'XDir','normal'); 
title('Delay-Time domain noise free Received Symbol Heat Map');
colormap(jet); colorbar;

figure;
imagesc(0:MN-1, 0:MN-1, abs(H_eq).^2);
set(gca,'XDir','normal'); 
title('DD domain Equivalent Channel Magnitude Heat Map');
colormap(jet); colorbar;

figure; bar3(abs(Y_DD).^2); title('DD domain pilot response Y\_DD (abs)'); xlabel('doppler'); ylabel('delay');
figure;
imagesc(0:M-1, 0:N-1, abs(Y_DD).^2);
set(gca,'XDir','normal'); 
title('DD domain pilot response Heat Map');
colormap(jet); colorbar;

figure; bar3(abs(Y_DD_noise_free).^2); title('DD domain noise free pilot response Y\_DD (abs)'); xlabel('doppler'); ylabel('delay');
figure;
imagesc(0:M-1, 0:N-1, abs(Y_DD_noise_free).^2);
set(gca,'XDir','normal'); 
title('DD domain noise free pilot response Heat Map');
colormap(jet); colorbar;
%{
figure;
imagesc(0:N-1, 0:N-1, real(F_N));
set(gca,'XDir','normal'); 
title('F\_N');
colormap(jet); colorbar;

figure;
imagesc(0:MN-1, 0:MN-1, real(F_NM));
set(gca,'XDir','normal'); 
title('F\_NM');
colormap(jet); colorbar;

%}