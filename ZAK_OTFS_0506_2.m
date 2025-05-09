% Combined DZT-based OTFS and 2-step OTFS Transceiver with DD-domain Pulse Shaping & explicit OTFS channel matrix
% All comments in English

clear; clc; close all;

%% ------------------ System Parameters ------------------
NUM_ITR   = 1;               % Number of Monte Carlo iterations for BER
B         = 0.96e6;          % Pulse bandwidth [Hz]
T         = 1.6e-3;          % Pulse duration [s]
nup       = 15e3;            % OTFS parameter nu_p [Hz]
taup      = 1/nup;           % Delay parameter tau_p [s]
M         = round(B * taup); % Number of delay bins
N         = round(T * nup);  % Number of Doppler bins
BT        = M * N;           % TF product (grid size) (MN)
symbvec   = (1/sqrt(2))*[1+1j,1-1j,-1+1j,-1-1j];  % 4-QAM constellation

%% ------------- Veh-A Channel Parameters -------------
numPaths   = 5;
delayP     = [0, 1, 2, 4, 7]/B;     % Path delays [s]
dopplerP   = [1, -2, -3, 3, 4]/T;   % Path Doppler [Hz]
powProfile = [0,-1,-9,-10,-13];    % Power profile [dB]
linpp      = 10.^(powProfile/10);
linpow     = linpp/sum(linpp);       % Normalized linear power

Fs = B;                  % Sampling rate = bandwidth
L  = M * N;              % Total time-domain samples
t  = (0:L-1)'/Fs;        % Time vector

%% ------------- Pre-generate symbols -------------
Zx_all = cell(NUM_ITR,1);
for itr = 1:NUM_ITR
    idx = randi([1,4], M*N, 1);
    Zx_all{itr} = reshape(symbvec(idx), M, N);
end

%% ----------------- SNR Sweep -----------------
SNRdB   = 10:15;
SNRlin  = 10.^(SNRdB/10);
BER_DZT = zeros(size(SNRdB));
BER_OTFS= zeros(size(SNRdB));

for si = 1:length(SNRdB)
    snr_lin = SNRlin(si);
    bitErr_D = 0; totalBits = 0;
    bitErr_O = 0;

    for itr = 1:NUM_ITR
        %% ---- 1) Generate multipath channel ----
        hgainP    = sqrt(linpow) .* exp(1j*2*pi*rand(1,numPaths));
        hgainmodP = hgainP .* exp(-1j*2*pi*delayP.*dopplerP);

        %% ---- 2) Build sparse DD-domain channel impulse response hdd ----
        hdd = zeros(2*M-1, 2*N-1);
        for p = 1:numPaths
            dp = round(delayP(p) * B);
            dq = round(dopplerP(p) * T);
            r  = M + dp;
            c  = N + dq;
            if r>=1 && r<=2*M-1 && c>=1 && c<=2*N-1
                hdd(r,c) = hdd(r,c) + hgainmodP(p);
            end
        end
        % 3D bar plot of DD-domain channel
        if itr ==1 && si ==1
        figure; bar3(abs(hdd));
        title('DD-domain channel impulse response |h_{dd}|');
        xlabel('Doppler index'); ylabel('Delay index'); zlabel('Magnitude');
        end

        %% ---- 3) Build explicit OTFS channel matrix in DD domain ----
        H_otfs = zeros(BT, BT);
        for kprime = 0:M-1
            for lprime = 0:N-1
                rx_idx = lprime*M + kprime + 1;
                for k = 0:M-1
                    for l = 0:N-1
                        tx_idx = l*M + k + 1;
                        tmp = 0;
                        for n = -1:1
                            for m = -1:1
                                ind1 = 2*M + kprime - k - n*M;
                                ind2 = 2*N + lprime - l - m*N;
                                if ind1>=1 && ind1<=2*M-1 && ind2>=1 && ind2<=2*N-1
                                    hval = hdd(ind1, ind2);
                                    tmp = tmp + hval * ...
                                        exp(1j*2*pi*(lprime - l - m*N)*(k + n*M)/BT) * ...
                                        exp(1j*2*pi*n*l/N);
                                end
                            end
                        end
                        H_otfs(rx_idx, tx_idx) = tmp;
                    end
                end
            end
        end

        %% ---- 4) Compute DD-domain MMSE equalizer ----
        MMSE_otfs = inv(H_otfs' * H_otfs + eye(BT)*(1/snr_lin)) * H_otfs';

        %% ---- Transmit: DZT-OTFS ----
        Zx   = Zx_all{itr};
        X_dzt = sqrt(N) * ifft(Zx, N, 2);
        x_dzt = X_dzt(:);

        %% ---- Transmit: 2-step OTFS ----
        Zx2    = Zx;
        X_tf   = sqrt(M*N) * ifft2(Zx2);
        S_ofdm = sqrt(M) * ifft(X_tf, M, 1);
        x_ofdm = S_ofdm(:);

        %% ---- Serial transmit frame plots ----
        figure;
        subplot(2,1,1);
        plot(real(x_dzt)); hold on; plot(imag(x_dzt));
        title('Time-domain transmit frame: DZT-OTFS');
        xlabel('Sample index'); legend('Real','Imag');
        subplot(2,1,2);
        plot(real(x_ofdm)); hold on; plot(imag(x_ofdm));
        title('Time-domain transmit frame: 2-step OTFS');
        xlabel('Sample index'); legend('Real','Imag');

        %% ---- Channel propagation & noise for DZT ----
        y_ch = zeros(L,1);
        for p = 1:numPaths
            di    = round(delayP(p)*Fs);
            xt_sh = [zeros(di,1); x_dzt(1:end-di)] .* exp(1j*2*pi*dopplerP(p)*t);
            y_ch  = y_ch + hgainmodP(p) * xt_sh;
        end
        noise = sqrt(1/(2*snr_lin))*(randn(L,1)+1j*randn(L,1));
        y_dzt = y_ch + noise;

        %% ---- Receiver: DZT equalization & DD-domain reconstruction ----
        Ydt      = reshape(y_dzt, M, N);
        Ydd      = fft(Ydt, N, 2) / sqrt(N);
        ydd_vec  = Ydd(:);
        xhat_dzt = MMSE_otfs * ydd_vec;
        Zhat_dzt = reshape(xhat_dzt, M, N);
        % 3D bar of equalized DD data (DZT)
        if itr ==1 && si ==1
        figure; bar3(abs(Zhat_dzt'));
        title('Equalized DD-domain data: DZT-OTFS');
        xlabel('Delay index'); ylabel('Doppler index'); zlabel('Magnitude');
        end
        %% ---- Receiver: 2-step OTFS equalization & DD-domain reconstruction ----
        y_ch2 = zeros(L,1);
        for p = 1:numPaths
            di     = round(delayP(p)*Fs);
            xt_sh2 = [zeros(di, 1); x_ofdm(1:end-di)] .* exp(1j*2*pi*dopplerP(p)*t);
            y_ch2  = y_ch2 + hgainmodP(p) * xt_sh2;
        end
        y_ofdm  = y_ch2 + noise;

        Y_tf    = reshape(y_ofdm, M, N);
        Y_tf    = fft(Y_tf, M, 1) / sqrt(M);
        Ydd2    = (1/sqrt(M*N)) * fft2(Y_tf);
        ydd2_vec= Ydd2(:);
        xhat_ofdm = MMSE_otfs * ydd2_vec;
        Zhat_ofdm = reshape(xhat_ofdm, M, N);
        % 3D bar of equalized DD data (2-step)
        if itr ==1 && si ==1
        figure; bar3(abs(Zhat_ofdm'));
        title('Equalized DD-domain data: 2-step OTFS');
        xlabel('Delay index'); ylabel('Doppler index'); zlabel('Magnitude');
        end
        %% ---- BER counting ----
        tx_vec      = Zx(:);
        [~, tx_id]  = min(abs(tx_vec*ones(1,4) - symbvec), [], 2);
        rx_vec_dzt  = Zhat_dzt(:);
        [~, rx_id]  = min(abs(rx_vec_dzt*ones(1,4) - symbvec), [], 2);
        rx_vec_ofdm = Zhat_ofdm(:);
        [~, rx2]    = min(abs(rx_vec_ofdm*ones(1,4) - symbvec), [], 2);

        bitErr_D = bitErr_D + sum(rx_id~=tx_id)*2;
        bitErr_O = bitErr_O + sum(rx2~=tx_id)*2;
        totalBits = totalBits + M*N*2;
    end

    BER_DZT(si)  = bitErr_D/totalBits;
    BER_OTFS(si) = bitErr_O/totalBits;
end

%% ----------------- Plot BER Curves -----------------
figure;
semilogy(SNRdB, BER_DZT, '-o','LineWidth',1.5); hold on;
semilogy(SNRdB, BER_OTFS,'-s','LineWidth',1.5);
grid on;
xlabel('SNR (dB)'); ylabel('BER');
legend('DZT-OTFS','2-step OTFS','Location','southwest');
title('BER Comparison of DZT vs 2-step OTFS with MMSE Equalization');
