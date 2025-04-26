clear all;
% =========================================================================
% Book title: OTFS Modulation: Theory and Applications
% Author: Saif Khan Mohammed, Ronny Hadani and Ananthanarayanan
% Chockalingam
% Editor: Wiley-IEEE Press
% Date: Oct. 2024
% Description: This Matlab program simulates the uncoded 4-QAM bit error
% rate (BER) for ZAK-OTFS (model-free) for the imperfect channel knowledge scenario.
% A separate Zak-OTFS frame is used for channel estimation. In the
% estimation frame a pilot (a single DD pulse at DD location (k_p, l_p))
% is transmitted and the channel response to this pulse is used to estimate
% h_{eff}[k,l] as discussed in equation (2.86) of the book.
% 4-QAM information symbols are then transmitted in a different Zak-OTFS frame
% whose response is then equalized using the estimate h_{eff}[k,l] in order
% to decode the information symbols.
% Veh-A channel is modeled. The program inputs are the SNR value for which
% the BER is to be simulated, bandwidth and time duration of Zak-OTFS pulsone,
% Zak-OTFS modulation parameter \nu_p, the maximum path Doppler shift \nu_{max}
% and the number of averaging iterations. DD domain sinc pulse shaping,
% i.e., w_rx(\tau, \nu) = w_{tx}(\tau, \nu) = \sqrt{B T} sinc(B \tau)
% sinc(T \nu) as given by equation (2.80) in the book.
% The averaged BER for ZAK-OTFS is saved in the variables
% berOTFS. This program can be used to simulate Zak-OTFS uncoded 4-QAM BER reported
% in Fig. 2.35, Fig. 3.16 and Fig. 3.17 of the book.
% =========================================================================

%----- Program INPUTs---------------------
NUM_ITR = 1000; % No. of iterations for BER averaging
B = 0.96e6; % Pulsone bandwidth (in Hz)
T = 1.6e-3; % Pulsone time (in seconds)
nup = 15e3; % Zak-OTFS modulation parameter \nu_p (in Hz)
taup = 1/nup;
SNRdB = 16; % SNR in dB
SNR_lin = 10^(SNRdB/10); % linear SNR
numax = 2000 ; % maximum Doppler shift (in Hz) induced by any path

%------- Channel parameters (Veh-A channel model)--------
numPaths = 6; % Six paths
delayP =  [ 0, 310, 710, 1090,  1730, 2510 ] * 1e-9 ; % Delay shifts induced by the 6 paths (in seconds)
powProfile = [0 ,-1, -9, -10, -15, -20]; % Average power gain of each path relative to the first path (in dB)
linpowProfile = (10.^(powProfile(1:numPaths)/10)) / (sum(10.^(powProfile(1:numPaths)/10))) ;



M = round(B*taup);
N = round(T*nup);
BT = round(B * T) ;
hgainP = zeros(numPaths,1);
symbvec = (1/sqrt(2))* [ (1+1i) , (1 - 1i), (-1 +1i), (-1-1i)] ; %vector of 4-QAM symbols

berOTFS=0;


for itr=1:1:NUM_ITR  


    for ijk=1:1:numPaths
    dopplerP(ijk) = numax*cos(2*pi*rand); % (0.5 - rand);
    end
    
        
        
    for pindex=1:1:numPaths
        %hgainP(pindex) = sqrt(linpowProfile(pindex)/2)*(randn + 1i*randn);
        hgainP(pindex) = sqrt(linpowProfile(pindex))*exp(1i*2*pi*rand);
    end

        
    
    btau = B*delayP;
    tdopplerP = T*dopplerP;
    


     % OTFS
    % generate h_dd[.,.]
    if (1)


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



    % Separate OTFS frames are used for channel estimation and data
    % transmission
    % Zak-OTFS Frame 1: A DD pulse is transmitted at (k_p, l_p) and its
    % response is used to estimate h_{eff}[k,l] (infinite pilot power is
    % assumed, i.e., there is no AWGN in the pilot frame)
    kp = round(M/2);
    lp = round(N/2);
    % variable hdd_est is the estimated \widehat{h}_{eff}[k,l]
    for kprime=0:1:(M-1)
        for lprime=0:1:(N-1)
           hdd_est(kp + kprime - kp + 1, lp + lprime -lp + 1) =  H_otfs(lprime*M + kprime + 1, lp*M + kp + 1) * exp(-1i*2*pi*(lprime - lp)*kp/BT);
        end
    end

    % H_otfs_est is the estimated effective DD domain channel matrix based
    % on the estimated channel filter hdd_est (constructed based on
    % equation (2.95) in book)

    H_otfs_est = zeros(BT, round(BT));
    for k=0:1:(M-1)
        for l=0:1:(N-1)
            for kprime=0:1:(M-1)
               for lprime=0:1:(N-1)
                  for n=-1:1:1
                      for m=-1:1:1
                          k_index = kp + kprime - k - n*M + 1 ;
                          l_index = lp + lprime - l -m*N + 1 ;
                          if ((k_index > 0) && (k_index <= M) && (l_index > 0) && (l_index <= N))
                              % this equation is based on (2.95) in the book
                             H_otfs_est(lprime*M + kprime + 1, l*M + k + 1) = hdd_est(k_index,l_index) *  exp(1i*2*pi*(lprime - l - m*N)*(k + n*M)/BT) * exp(1i*2*pi*n*l/N) ;
                          end
                      end
                  end
               end
            end
        end
    end



     
    %Zak-OTFS MMSE matrix for equalization at receiver
    MMSE_otfs = inv(H_otfs_est'*H_otfs_est + eye(BT,BT)*(1/SNR_lin)) * H_otfs_est' ;
    end

    % Information symbols are transmitted in a different Zak-OTFS frame
    xmsg = randi([1 4 ], M*N, 1);
    xinft = zeros(M*N,1);
    for u=1:1:M*N
     xinft(u) = symbvec(xmsg(u)) ;
    end
    
    % H_otfs * xinft is the vector of received DD domain samples
    y_otfs = H_otfs * xinft  +  sqrt(1/(2*SNR_lin))*(randn(BT,1) + 1i*randn(BT,1)) ;
    
    % MMSE equalization of the received DD samples using the estimated
    % channel
    xrx_otfs = MMSE_otfs * y_otfs ;
    
    % Accumulation of bit errors
    berOTFS = berOTFS + length(find(sign(real(xinft)) ~= sign(real(xrx_otfs)))) + length(find(sign(imag(xinft)) ~= sign(imag(xrx_otfs))));
    
    
      
    
    disp('------------------------');
    disp(itr);
    disp('BER Zak-OTFS');
    disp(berOTFS/(itr*BT*2));


end


berOTFS = berOTFS/(NUM_ITR*BT*2);


    



