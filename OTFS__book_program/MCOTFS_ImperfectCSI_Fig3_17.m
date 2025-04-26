clear all;

% =========================================================================
% Book title: OTFS Modulation: Theory and Applications
% Author: Saif Khan Mohammed, Ronny Hadani and Ananthanarayanan
% Chockalingam
% Editor: Wiley-IEEE Press
% Date: Oct. 2024
% Description: This Matlab program simulates the uncoded 4-QAM bit error
% rate (BER) for Multicarrier (MC)-OTFS (model-free) for the imperfect channel knowledge scenario.
% A separate Zak-OTFS frame is used for channel estimation. In the
% estimation frame a pilot (a single DD pulse at DD location (k_p, l_p))
% is transmitted and the channel response to this pulse is used to estimate
% the DD domain I/O relation.
% 4-QAM information symbols are then transmitted in a different MC-OTFS frame
% whose response is then equalized using the channel estimate in order
% to decode the information symbols.
% Veh-A channel is modeled. The program inputs are the SNR value for which
% the BER is to be simulated, bandwidth and time duration of MC-OTFS pulsone,
% MC-OTFS modulation parameter \nu_p, the maximum path Doppler shift \nu_{max}
% and the number of averaging iterations. Rectangular time-frequency
% window W_{tx}(t,f) is used as in equation (3.23) of the book.
% The averaged BER for MC-OTFS is saved in the variables
% berOTFS. This program can be used to simulate MC-OTFS uncoded 4-QAM BER reported
% in Fig. 3.16 and Fig. 3.17 of the book.
% =========================================================================

%----- Program INPUTs---------------------
NUM_ITR = 10; % No. of iterations for BER averaging
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
L = numPaths;



taunorm = delayP/taup;



for itr=1:1:NUM_ITR
    
    for ijk=1:1:numPaths
    dopplerP(ijk) = numax*cos(2*pi*rand);
    end
    
    nunorm = dopplerP/nup;
    
    for pindex=1:1:numPaths
        
        hgainP(pindex) = sqrt(linpowProfile(pindex))*exp(1i*2*pi*rand);
    end    
    
    % Calculate the DD domain I/O relation function h_w[k,l;k',l'] as
    % defined in equation (3.90) of the book
    % Variable H models the I/O relation function h_w in the book
    % h_w[k,l ; k',l'] = \sum_i H[i; l+k*M +1, l'+k'*M + 1]
    H = zeros(numPaths, M*N, M*N);
    % Sinc table 1
    for i=1:1:L
        for mpmm=-(M-1):1:(M-1)
            sinctab1(i,M+mpmm) = sinc((1 - taunorm(i)) * ( nunorm(i) - mpmm) ) ;
        end
    end
    
    %exptable 1
    
    
    for i=1:1:L
        for mpmm=-(M-1):1:(M-1)
            tmmp = nunorm(i) - (mpmm) ;
            extab1(i,M+mpmm) =  exp(1i*pi*taunorm(i)*tmmp) ;
        end
    end
     
    % Sinc Table 2
    for i=1:1:L
        for mpmm=-(M-1):1:(M-1)
            tmmp = nunorm(i) - (mpmm) ;
            sinctab2A(i,M+mpmm) = (1 - taunorm(i))*exp(1i*pi*tmmp) * sinc((1 - taunorm(i))*tmmp)  ;
        end
    end
    
        for i=1:1:L
        for mpmm=-(M-1):1:(M-1)
            tmmp = nunorm(i) - (mpmm) ;
            sinctab2B(i,M+mpmm) =     taunorm(i)*sinc(taunorm(i)*tmmp) ;
        end
        end
    
        
    
    
    for m_id=0:1:(M-1)
        for l_id=0:1:(M-1)
            exptab22(m_id+1,l_id+1) = exp(-1i*2*pi*m_id*l_id/M) ;
        end
        
        for i=1:1:L
            exptab33(m_id+1,i) = exp(1i* 2*pi*(  - 0.5*( - m_id)*(1 - taunorm(i))   )   ) ;
        end
    end
    
    
                    for ijk=1:1:L
                        for m_id=0:1:(M-1)
                            tmtab1(ijk,m_id+1) = exp(-1i*2*pi*m_id* taunorm(ijk));
                            for l=0:1:M-1
                            tmtab2(l+1,m_id+1) = exp(-1i*2*pi*m_id*l/M) ;
                            end
                        end
                    end
    
    
    for kp=0:1:(N-1)
        
        
        for lp=0:1:(M-1)
            
            kplp = lp + kp*M + 1;
            
            
 
                
                for l=0:1:(M-1)
                    
  
                    tl =0;
                    
                    
                    
                    for i=1:1:L
                        
                        % term for n'=0
                        htmp =0;
                        
                            mp = 0:1:(M-1) ;
                            thetamp = exp(1i*2*pi*mp*(lp/M - 0.5*taunorm(i) -0.5)) ;
                            
                            
                            for m=0:1:(M-1)
                            
                            
                            sincvec1 = sinctab1(i,M + mp - m);
                            
                            htmp = htmp + exptab22(m+1,l+1)*exptab33(m+1,i)*sum( thetamp .* sincvec1 );
                            end
                        hzero = htmp*(1 - taunorm(i))*exp(1i*pi*(1 - taunorm(i))*nunorm(i)) ;
                        
                        % term for other n'=1 .. N-1
                        h1A =0;
                        h1B =0;
                        mp = 0:1:(M-1) ;
                        expvec11 = exp(1i*2*pi*mp*lp/M) ;
                        
                        
                        for m=0:1:(M-1)
                            
                            tm = tmtab1(i,m+1)*tmtab2(l+1,m+1);
                           
                            exvec1 = extab1(i,M+mp - m);
                            sincvec2A = sinctab2A(i,M + mp - m);
                            sincvec2B = sinctab2B(i,M + mp - m);
                            h1A = h1A + tm*sum(expvec11 .* exvec1 .* sincvec2A  ) ;
                            h1B = h1B + tm*sum(expvec11 .* exvec1 .* sincvec2B  ) ;
                        end
                    
                    for k=0:1:(N-1)
                        kl = l + k*M +1 ;
                        ab = nunorm(i) + (k - kp)/N;
                        h1 = (h1A + exp(-1i*2*pi*k/N)*h1B) * exp(-1i*2*pi*nunorm(i)*taunorm(i))*exp(1i*2*pi*ab)*(N-1)*exp(1i*pi*(N-2)*ab)*sinc((N-1)*ab)/sinc(ab) ;
                        
                        H(i,kplp, kl) = hgainP(i)*(hzero + h1)/(M*N) ;
                        
                    end % for k=0:1:(N-1)
                    
                    end % for i=1:1:L   
                    
                    
                    
                end
           
        end
    end
   
    H_otfs = zeros(BT,BT);
    for ijk=1:1:numPaths
    H_otfs = H_otfs + reshape(H(ijk,:,:), BT, BT);
    end
    clear H;
    
    % Separate MC-OTFS frames are used for channel estimation and data
    % transmission
    % MC-OTFS Frame 1: A DD pulse is transmitted at (k_p, l_p) and its
    % response is used to estimate h_w[k,l ; k',l'] (infinite pilot power is
    % assumed, i.e., there is no AWGN in the pilot frame)
    kp = N/2;
    lp = M/2;
    hdd_est = zeros(N,M);

    for kprime=0:1:(N-1)
        for lprime=0:1:(M-1)
            hdd_est(kprime +1 , lprime+1) = H_otfs(lprime + kprime*M + 1, lp + kp*M + 1);
        end
    end

    H_otfs_est = zeros(M*N, M*N);
    for kprime=0:1:(N-1)
        for lprime=0:1:(M-1)
            for k=0:1:(N-1)
                for l=0:1:(M-1)
                    kpt =  mod(kp + kprime -k,N) + 1;
                    lpt =  mod(lp + lprime -l,M) + 1;
                    if ((kpt > 0) && (kpt <= N) && (lpt > 0) && (lpt <= M))
                       H_otfs_est(lprime + kprime*M + 1, l+k*M + 1) = hdd_est(kpt, lpt );
                    end
                end
            end
        end
    end


    %MC-OTFS MMSE matrix for equalization at receiver
    MMSE_mcotfs = inv(H_otfs_est'*H_otfs_est + eye(BT,BT)*(1/SNR_lin)) * H_otfs_est' ;
    
    % Information symbols are transmitted in a different MC-OTFS frame
    xmsg = randi([1 4 ], M*N, 1);
    xinft = zeros(M*N,1);
    for u=1:1:M*N
     xinft(u) = symbvec(xmsg(u)) ;
    end
    
    % H_otfs * xinft is the vector of received DD domain samples
    y_otfs = H_otfs * xinft  +  sqrt(1/(2*SNR_lin))*(randn(BT,1) + 1i*randn(BT,1)) ;
    
    % MMSE equalization of the received DD samples using the estimated
    % channel
    xrx_otfs = MMSE_mcotfs * y_otfs ;
    
    % Accumulation of bit errors
    berOTFS = berOTFS + length(find(sign(real(xinft)) ~= sign(real(xrx_otfs)))) + length(find(sign(imag(xinft)) ~= sign(imag(xrx_otfs))));
    
    disp(berOTFS/(itr*2*BT));

end

berOTFS = berOTFS/(NUM_ITR*2*BT);
    




  


