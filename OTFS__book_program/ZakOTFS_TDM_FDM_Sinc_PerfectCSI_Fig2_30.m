clear all;
% =========================================================================
% Book title: OTFS Modulation: Theory and Applications
% Author: Saif Khan Mohammed, Ronny Hadani and Ananthanarayanan
% Chockalingam
% Editor: Wiley-IEEE Press
% Date: Oct. 2024
% Description: This Matlab program simulates the uncoded 4-QAM bit error
% rate (BER) for TDM, FDM and ZAK-OTFS for the perfect channel knowledge scenario.
% Veh-A channel is modeled. The program inputs are the SNR value for which
% the BER is to be computed, bandwidth and time duration of Zak-OTFS pulsone,
% Zak-OTFS modulation paraneter \nu_p and number of averaging iterations.
% The averaged BER for TDM, FDM and ZAK-OTFS is saved in the variables
% berTDM, berFDM and berOTFS respectively. This program simulates a single
% SNR point for the curves in Fig. 2.30 of the book.
% =========================================================================

%----- Program INPUTs---------------------
NUM_ITR = 1000; % No. of iterations for BER averaging
B = 0.96e6; % Pulsone bandwidth in Hz
T = 1.6e-3; % Pulsone time in seconds
nup = 15e3; % Zak-OTFS modulation parameter \nu_p in Hz
taup = 1/nup;
SNRdB = 15; % SNR in dB
SNR_lin = 10^(SNRdB/10); % linear SNR

%------- Channel parameters (Veh-A channel model)--------
numPaths = 6; % Six paths
delayP =  [ 0, 310, 710, 1090,  1730, 2510 ] * 1e-9 ; % Delay shifts induced by the 6 paths in seconds
powProfile = [0 ,-1, -9, -10, -15, -20]; % Average power gain of each path relative to the first path (in dB)
linpowProfile = (10.^(powProfile(1:numPaths)/10)) / (sum(10.^(powProfile(1:numPaths)/10))) ;
numax = 815 ; % maximum Doppler shift induced by any path in Hz


M = round(B*taup);
N = round(T*nup);
BT = round(B * T) ;
hgainP = zeros(numPaths,1);
symbvec = (1/sqrt(2))* [ (1+1i) , (1 - 1i), (-1 +1i), (-1-1i)] ; %vector of 4-QAM symbols
berTDM  = 0;
berOTFS=0;
berFDM = 0;



for itr=1:1:NUM_ITR  

    for ijk=1:1:numPaths
    dopplerP(ijk) = numax *cos(2*pi*rand);
    end
        
    for pindex=1:1:numPaths
        hgainP(pindex) = sqrt(linpowProfile(pindex))*exp(1i*2*pi*rand);
    end
    
    
    btau = B*delayP;
    tdopplerP = T*dopplerP;

    
    
    
    %H_mat: I/O relation matrix for TDM
    H_mat = zeros(round(4+B*(T + delayP(numPaths))), BT);
    
    for kprime=1:1:round(4+B*(T + delayP(numPaths)))
        for pindex=1:1:numPaths
        for k=max(1,kprime - 20 - round(btau(pindex))):1:min(BT, kprime + 20 - round(btau(pindex))) 
            H_mat(kprime, k) = H_mat(kprime, k) +  hgainP(pindex)* exp(1i*2*pi*dopplerP(pindex)*(k-1)/B) * (1 - abs(dopplerP(pindex))/B)  *  exp(1i*pi*dopplerP(pindex)*(kprime - k - btau(pindex))/B) *  sinc( (1 - abs(dopplerP(pindex))/B)*(kprime - k - btau(pindex))) ;
        end
        end
    end
    MMSE_mat = inv(H_mat'*H_mat + eye(BT,BT)*(1/SNR_lin)) * H_mat' ;
    % End of generation of TDM I/O matrix
    
    
    
    % Zak-OTFS
    
    for ijk=1:1:numPaths
    hgainmodP(ijk) =  hgainP(ijk)*exp(-1i*2*pi*delayP(ijk)*dopplerP(ijk));
    end
    % Variable hdd(.,.) is the effective DD domain channel h_{eff}[k,l] with
    % sinc pulse shaping

    hdd = zeros(4*M-1, 4*N-1);
    for ind1=1:1:(4*M-1)
        for ind2 =1:1:(4*N-1)
            p = ind1 - 2*M ;
            q = ind2 - 2*N ;
            tmp_sum = 0;
            
            for ijk=1:1:numPaths
                f12 = @(x) exp(1i*pi*(p + btau(ijk))*x/BT) .* sinc(q - x)  .* sinc(x - tdopplerP(ijk)) .* (1 - abs(x/BT)) .* sinc((p - btau(ijk)) * (1 - abs(x/BT)) )  ;
            
                tmp_sum = tmp_sum + hgainmodP(ijk)*integral(f12, max(max(-BT,tdopplerP(ijk) - 20), q - 20), min(min(BT,tdopplerP(ijk) + 20), q+20));
            end
            hdd(ind1, ind2) = tmp_sum ;
        end
    end
    
    % H_otfs is the effective MN X MN DD domain channel matrix, see the
    % mattrix H_{dd} in equations (2.95) and (2.96) of the book
    % equation (2.96) describes how the elements of H_{d} matrix depend on
    % the effective DD domain channel fikter h_{eff}[k,l]    

    H_otfs = zeros(BT, round(BT)); 
    for kprime=0:1:(M-1)
        for lprime=0:1:(N-1)
            for k=0:1:(M-1)
                for l=0:1:(N-1)
                    for n=-1:1:1
                        for m=-1:1:1
                    
                     H_otfs(lprime*M + kprime + 1, l*M + k + 1) =  H_otfs(lprime*M + kprime + 1, l*M + k + 1) + hdd(2*M + kprime - k - n*M, 2*N+lprime - l - m*N) * exp(1i*2*pi*(lprime - l - m*N)*(k + n*M)/BT) * exp(1i*2*pi*n*l/N);
                        
                        end
                    end
                end
            end
        end
    end
    
    % compute MMSE matrix for equalization at receiver
    MMSE_otfs = inv(H_otfs'*H_otfs + eye(BT,BT)*(1/SNR_lin)) * H_otfs' ;
   

    % FDM I/O matrix generation
    kprime_min = 0 - round(max(abs(dopplerP(1:numPaths)))*T) - 4;
    kprime_max = (BT-1) + round(max(abs(dopplerP(1:numPaths)))*T) + 4;
    
    hfd = zeros(2*(kprime_max - kprime_min)+1,BT);
    for ijk=1:1:numPaths
    for n=(round(tdopplerP(ijk)) - 20):1:(round(tdopplerP(ijk)) + 20)   
        for k=0:1:(BT-1)
            hfd(n + (kprime_max - kprime_min) + 1, k+1) =  hfd(n + (kprime_max - kprime_min) + 1, k+1) +   hgainP(ijk)* (1 - delayP(ijk)/T) * exp(-1i*2*pi*k*delayP(ijk)/T) * exp(-1i*pi*(n + tdopplerP(ijk))*delayP(ijk)/T) * sinc((n - tdopplerP(ijk)) * (1 - delayP(ijk)/T)) ; ;
        end
    end
    end

    H_fdm = zeros(kprime_max - kprime_min + 1, BT);
    
    for kprime=kprime_min:1:kprime_max
        for k=0:1:(BT-1)
            H_fdm(kprime+1-kprime_min, k+1) = hfd(kprime - k + (kprime_max - kprime_min) + 1, k+1) ;
        end 
    end
    MMSE_fdm = inv(H_fdm'*H_fdm + eye(BT,BT)*(1/SNR_lin)) * H_fdm' ;
    % End of FDM I/O matrix generation
    
    
  
    % Transmit TDM frame, equalize it and compute accumulate bit errors
    xmsg = randi([1 4 ], M*N, 1);
    xinf = zeros(M*N,1);
    for u=1:1:M*N
     xinf(u) = symbvec(xmsg(u)) ; % vector of transmitted symbols
    end
    % received vector of samples in TDM
    yrx = H_mat*xinf + sqrt(1/(2*SNR_lin))*(randn(round(4+B*(T + delayP(numPaths))),1) + 1i*randn(round(4+B*(T + delayP(numPaths))),1)) ;
    % Equaization of receiuved TDM samples
    xrx_hat = MMSE_mat * yrx ;
    % Accumulate of bit errors
    berTDM = berTDM + length(find(sign(real(xinf)) ~= sign(real(xrx_hat)))) + length(find(sign(imag(xinf)) ~= sign(imag(xrx_hat))));
    % End of simulation of TDM frame   
       
    % Simulating a Zak-OTFS frame transmission, reception and equalization
    xmsg = randi([1 4 ], M*N, 1);
    xinft = zeros(M*N,1);
    for u=1:1:M*N
     xinft(u) = symbvec(xmsg(u)) ; % Vector of MN DD domain 4-QAM information symbols
    end
    
    % H_otfs * xinft is the received DD domain samples
    y_otfs = H_otfs * xinft  +  sqrt(1/(2*SNR_lin))*(randn(BT,1) + 1i*randn(BT,1)) ;
    % MMSE receiver equalization    
    xrx_otfs = MMSE_otfs * y_otfs ;
    % Accumulation of bit errors
    berOTFS = berOTFS + length(find(sign(real(xinft)) ~= sign(real(xrx_otfs)))) + length(find(sign(imag(xinft)) ~= sign(imag(xrx_otfs))));
    
    
    % FDM< frame generation, transmission, reception and equalization
    xmsg = randi([1 4 ], M*N, 1);
    xinf = zeros(M*N,1);
    for u=1:1:M*N
     xinf(u) = symbvec(xmsg(u)) ; % Frequency domain (FD) 4-QAM information symbols
    end
    
    % H_fdm*xinf is the received FD samples
    yfdm = H_fdm*xinf + sqrt(1/(2*SNR_lin))*(randn(kprime_max - kprime_min + 1,1) + 1i*randn(kprime_max - kprime_min + 1,1)) ;
    
    % MMSE equalization in FD
    xfdm_hat = MMSE_fdm * yfdm ;
    
    berFDM = berFDM + length(find(sign(real(xinf)) ~= sign(real(xfdm_hat)))) + length(find(sign(imag(xinf)) ~= sign(imag(xfdm_hat))));

    disp('------------------------');
    disp(itr);
    disp('TDM');
    disp(berTDM/(itr*BT * 2));
    disp('OTFS');
    disp(berOTFS/(itr*BT * 2));
    disp('FDM');
    disp(berFDM/(itr*BT * 2));
    

    
    
end % for itr=1:1:NUM_ITR
    
    
% Computation of average BER for TDM, FDM and Zak-OTFS (perfect CSI)
berTDM = berTDM/(NUM_ITR*BT * 2);
berOTFS = berOTFS/(NUM_ITR*BT*2);
berFDM = berFDM/(NUM_ITR*BT * 2);




