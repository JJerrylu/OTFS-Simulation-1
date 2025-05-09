clear; clc; close all;

%% ----- Program INPUTs （保留原始參數） -----
NUM_ITR   = 10;               % BER averaging 次數
B         = 0.96e6;            % 脈波帶寬 [Hz]
T         = 1.6e-3;            % 脈波時間 [s]
nup       = 15e3;              % Zak‑OTFS 多普勒參數 ν_p [Hz]
taup      = 1/nup;             % 延遲參數 τ_p
M         = round(B * taup);   % delay bins
N         = round(T * nup);    % doppler bins
BT        = round(B * T);      % = M * N （time–freq 點數）
symbvec   = (1/sqrt(2))*[1+1j,1-1j,-1+1j,-1-1j];  % QPSK

SNRdB_vec = 0:1:15;
SNR_lin   = 10.^(SNRdB_vec/10);
BER       = zeros(size(SNRdB_vec));

% ----- channel parameters -----
numPaths    = 6;
delayP      = [ 0, 310, 710, 1090, 1730, 2510 ] * 1e-9;          % [s]
powProfile  = [0, -1, -9, -10, -15, -20];                       % [dB]
linpp       = 10.^(powProfile/10);
linpowProfile = linpp / sum(linpp);
numax       = 815;                                              % [Hz]

Fs = B;                % 採樣頻率 = 脈波帶寬
L  = M * N;            % 總樣本數
t  = (0:L-1)'/Fs;      % 時域對應向量

% 產生多徑參數
dopplerP = numax * cos(2*pi*rand(1,numPaths));   % 隨機多普勒
hgainP   = sqrt(linpowProfile) .* exp(1j*2*pi*rand(1,numPaths));  % 路徑增益
hgainmodP= hgainP .* exp(-1j*2*pi*delayP.*dopplerP); % 考慮延遲-多普勒相位

%% ----- 預先產生所有迭代的 DD‑domain 資料 -----
Zx_all = cell(NUM_ITR,1);
for itr = 1:NUM_ITR
  msg = randi([1 4], M*N, 1);
  xdd = symbvec(msg);
  % reshape 成 M×N (delay×doppler)
  Zx_all{itr} = reshape(xdd, M, N);
end

%% ----- SNR 迴圈 -----
for si = 1:numel(SNRdB_vec)
  snr_lin = SNR_lin(si);
  bitErr = 0;
  totalBits = 0;

  for itr = 1:NUM_ITR
    Zx = Zx_all{itr};   % M×N 矩陣

    %% 1) TX: IDZT (DD → delay–time)
    x_t = zeros(M*N,1);
    for n = 0:(M*N-1)
      mu = floor(n/M);
      ell = mod(n, M);
      acc = 0;
      for k = 0:(N-1)
        acc = acc + Zx(ell+1, k+1) * exp(1j*2*pi * mu * (k/N));
      end
      x_t(n+1) = acc / sqrt(N);
    end

    %% 第一次畫圖：TX 端的 3D DD-domain & delay-time domain
    if itr==1 && si==1
      figure; bar3(abs(Zx.') ); 
      title('TX DD–domain |Z_x[k,\ell]|'); xlabel('Delay \ell'); ylabel('Doppler k');
      
      figure; 
      X_dt = reshape(x_t, M, N);   % M rows (= delay bins), N cols (= time segment)
      bar3(abs(X_dt.') ); 
      title('TX delay–time |x_t[n]|'); xlabel('Time seg. \mu'); ylabel('delay \ell');
    end
    
    y_ch = zeros(L,1);
    for i = 1:numPaths
        di = round(delayP(i)*Fs);            % 延遲對應的樣本點
        % 對 x_t 做延遲、頻移
        xt_shifted = [zeros(di,1); x_t(1:end-di)] .* exp(1j*2*pi*dopplerP(i)*t);
        y_ch = y_ch + hgainmodP(i) * xt_shifted;
    end

    
    %% 2) AWGN 通道
    noise = sqrt(1/(2*snr_lin))*(randn(size(x_t)) + 1j*randn(size(x_t)));
    y_t   = y_ch + noise;

    %% 3) RX S/P → 重組成 delay–time domain 格點 M×N
    Y_dt = reshape(y_t, M, N);

    if itr==1 && si==1
      figure; bar3(abs(Y_dt.') );
      title('RX delay–time |Y_{dt}[\ell,\mu]|'); xlabel('Time seg. \mu'); ylabel('delay \ell');
    end

    %% 4) DZT (delay–time → DD)
    Y_DD = zeros(M,N);
    for ell = 0:(M-1)
      for k = 0:(N-1)
        acc = 0;
        for mu = 0:(N-1)
          acc = acc + Y_dt(ell+1, mu+1) * exp(-1j*2*pi * k * mu / N);
        end
        Y_DD(ell+1, k+1) = acc / sqrt(N);
      end
    end

    if itr==1 && si==1
      figure; bar3(abs(Y_DD.') );
      title('RX DD–domain |Y_{DD}[k,\ell]| after DZT'); xlabel('Delay \ell'); ylabel('Doppler k');
    end

    %% 5) Equalization (AWGN, perfect CSI → identity)
    Xhat = Y_DD;

    %% 6) Demapping & BER
    xhat_vec = Xhat(:);
    tx_vec   = reshape(Zx, [], 1);
    % 對每個符號做最接近映射
    [~, rx_idx] = min(abs(xhat_vec*ones(1,4) - ones(M*N,1)*symbvec), [], 2);
    [~, tx_idx] = min(abs(tx_vec*ones(1,4)   - ones(M*N,1)*symbvec), [], 2);
    % QPSK 每符號 2 bits
    bitErr   = bitErr   + sum(rx_idx~=tx_idx)*2;
    totalBits= totalBits+ M*N*2;
  end

  BER(si) = bitErr / totalBits;
end

%% 7) 繪製 SNR–BER 曲線
figure;
semilogy(SNRdB_vec, BER, '-o','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('DZT–QPSK over AWGN (100 itrs)');
