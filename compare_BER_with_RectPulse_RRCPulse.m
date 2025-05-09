% compare_ber_otfs.m
% This script compares the BER performance of 2-step OTFS (ISFFT+OFDM) vs DZT-OTFS
% under two pulse-shaping scenarios: RRC and ideal rectangular.

%% Parameters
M = 16;          % number of Doppler bins
N = 16;          % number of delay bins
numIter = 1000;  % Monte Carlo iterations per SNR point
snrVec = 0:2:20; % SNR range (dB)

%% 1) RRC Pulse Shaping
rolloff = 0.3; span = 6;      % RRC parameters
pulse_rrc = rcosdesign(rolloff, span, N, 'sqrt');
ber2_rrc = zeros(size(snrVec));
berD_rrc = zeros(size(snrVec));

for idx = 1:length(snrVec)
    errors2 = 0; errorsD = 0; totalBits = 0;
    for iter = 1:numIter
        %--- Transmitter ---
        bits = randi([0 1], M*N*2, 1);
        symbols = qpskMod(bits);
        ddGrid = reshape(symbols, M, N);
        tx2 = otfs2step(ddGrid, pulse_rrc);
        txD = otfsDZT(ddGrid, pulse_rrc);

        %--- Channel ---
        rx2 = awgn(tx2, snrVec(idx), 'measured');
        rxD = awgn(txD, snrVec(idx), 'measured');

        %--- Receiver ---
        ddRec2 = otfsRx2step(rx2, pulse_rrc, M, N);
        ddRecD = otfsRxDZT(rxD, pulse_rrc, M, N);

        %--- Detection ---
        bits2 = qpskDemod(ddRec2(:));
        bitsD = qpskDemod(ddRecD(:));
        errors2 = errors2 + sum(bits2 ~= bits);
        errorsD = errorsD + sum(bitsD ~= bits);
        totalBits = totalBits + length(bits);
    end
    ber2_rrc(idx) = errors2/totalBits;
    berD_rrc(idx) = errorsD/totalBits;
end

figure;
semilogy(snrVec, ber2_rrc, '-o', snrVec, berD_rrc, '-s');
grid on;
xlabel('SNR (dB)'); ylabel('BER');
legend('2-step OTFS','DZT-OTFS','Location','best');
title('BER under RRC Pulse Shaping');

%% 2) Rectangular Pulse Shaping
pulse_rect = ones(N,1);   % ideal rectangular pulse of length N
ber2_rect = zeros(size(snrVec));
berD_rect = zeros(size(snrVec));

for idx = 1:length(snrVec)
    errors2 = 0; errorsD = 0; totalBits = 0;
    for iter = 1:numIter
        bits = randi([0 1], M*N*2, 1);
        symbols = qpskMod(bits);
        ddGrid = reshape(symbols, M, N);
        tx2 = otfs2step(ddGrid, pulse_rect);
        txD = otfsDZT(ddGrid, pulse_rect);
        rx2 = awgn(tx2, snrVec(idx), 'measured');
        rxD = awgn(txD, snrVec(idx), 'measured');
        ddRec2 = otfsRx2step(rx2, pulse_rect, M, N);
        ddRecD = otfsRxDZT(rxD, pulse_rect, M, N);
        bits2 = qpskDemod(ddRec2(:));
        bitsD = qpskDemod(ddRecD(:));
        errors2 = errors2 + sum(bits2 ~= bits);
        errorsD = errorsD + sum(bitsD ~= bits);
        totalBits = totalBits + length(bits);
    end
    ber2_rect(idx) = errors2/totalBits;
    berD_rect(idx) = errorsD/totalBits;
end

figure;
semilogy(snrVec, ber2_rect, '-o', snrVec, berD_rect, '-s');
grid on;
xlabel('SNR (dB)'); ylabel('BER');
legend('2-step OTFS','DZT-OTFS','Location','best');
title('BER under Rectangular Pulse Shaping');

%% Supporting Functions
function txSig = otfs2step(ddGrid, pulse)
    [M, N] = size(ddGrid);
    % ISFFT
    tfGrid = ifft(ifft(ddGrid.').', M, 1) * sqrt(M*N);
    % OFDM modulator with pulse shaping
    txSym = [];
    for n = 1:N
        ofdmSym = ifft(tfGrid(:,n), M);
        txSym = [txSym; conv(ofdmSym, pulse, 'same')];
    end
    txSig = txSym;
end

function ddRec = otfsRx2step(rxSig, pulse, M, N)
    rxGrid = zeros(M, N);
    ptr = 1;
    for n = 1:N
        segment = rxSig(ptr:ptr+M-1);
        segment = conv(segment, pulse, 'same');
        rxGrid(:,n) = fft(segment, M);
        ptr = ptr + M;
    end
    ddRec = fft(fft(rxGrid).').'/sqrt(M*N);
end

function txSig = otfsDZT(ddGrid, pulse)
    [M, N] = size(ddGrid);
    z = reshape(ddGrid, M*N, 1);
    txSig = conv(z, pulse, 'same');
end

function ddRec = otfsRxDZT(rxSig, pulse, M, N)
    y = conv(rxSig, pulse, 'same');
    ddRec = reshape(y(1:M*N), M, N);
end

function sym = qpskMod(bits)
    b = reshape(bits, 2, []).';
    sym = (1/sqrt(2))*((1-2*b(:,1)) + 1j*(1-2*b(:,2)));
end

function bits = qpskDemod(sym)
    b = zeros(length(sym),2);
    b(:,1) = real(sym) < 0;
    b(:,2) = imag(sym) < 0;
    bits = reshape(b.', [], 1);
end
