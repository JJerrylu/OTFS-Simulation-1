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
% and the number of averaging iterations. DD domain RRC pulse shaping
% is used as given by equation (2.81) in the book.
% The averaged BER for ZAK-OTFS is saved in the variables
% berOTFS. This program can be used to simulate Zak-OTFS uncoded 4-QAM BER reported
% in Fig. 2.36, Fig. 3.16 and Fig. 3.17 of the book.
% =========================================================================


%---------- INPUT Parameters ------------------------------------
NUM_ITR = 10; % number of averaging iterations
B = 0.96e6; % Pulsone bandwidth in Hz
T = 1.6e-3; % Pulsone time duration in seconds
nup = 30e3; % Zak-OTFS modulation parameter \nu_p (in Hz)
taup = 1/nup;
beta_nu = 0.2; % RRC roll-off factor \beta_nu in book for pulse shape along Doppler domain
beta_tau = 0.1; % RRC roll-off factor \beta_nu in book for pulse shape along delay domain
SNRdB = 16; % Signal to noise ratio (in dB) 
SNR_lin = 10^(SNRdB/10);
numax = 2000 ; % maximum Doppler shift (in Hz) induced by any path

%----- Channel path parameters (Veh-A) --------------------------
numPaths = 6;
delayP =  B*[ 0, 310, 710, 1090,  1730, 2510 ] * 1e-9 ;
powProfile = [0 ,-1, -9, -10, -15, -20]; % Average power gain of each path relative to the first path (in dB)
linpowProfile = (10.^(powProfile(1:numPaths)/10)) / (sum(10.^(powProfile(1:numPaths)/10))) ;


symbvec = (1/sqrt(2))* [ (1+1i) , (1 - 1i), (-1 +1i), (-1-1i)] ; % vector of 4-QAM information symbols



% Zak-OTFS modulation parameters
M = round(B*taup);
N = round(T*nup);
BT = round(B * T) ;

hgainP = zeros(numPaths,1);
berOTFS=0;


constbetanu = (beta_nu*sqrt(1/2))* ( (1 + (2/pi))*sin(pi/(4*beta_nu))   +  (1 - (2/pi))*cos(pi/(4*beta_nu))) ;
abstolnu = 0.001;
constbetatau = (beta_tau*sqrt(1/2))* ( (1 + (2/pi))*sin(pi/(4*beta_tau))   +  (1 - (2/pi))*cos(pi/(4*beta_tau))) ;
abstoltau = 0.001;

for itr=1:1:NUM_ITR  


    for ijk=1:1:numPaths
    dopplerP(ijk) = numax*cos(2*pi*rand); % (0.5 - rand);
    end
    
    
        
        
    for pindex=1:1:numPaths
        hgainP(pindex) = sqrt(linpowProfile(pindex))*exp(1i*2*pi*rand);
    end

    dopplerP = dopplerP * T;
    

    % Variable hdd(.,.) is the effective DD domain channel h_{eff}[k,l] with
    % RRC pulse shaping

    hdd = zeros(4*M-1, 4*N-1);
    for ijk=1:1:numPaths
        if (dopplerP(ijk) > 0)
            idopplerP = ceil(dopplerP(ijk));
        else
            idopplerP = floor(dopplerP(ijk));
        end
        
       for ind1=max(1,2*M + ceil(delayP(ijk)) - 20 ):1:min((4*M-1),2*M + ceil(delayP(ijk)) + 20 ) 
           for ind2 =max(1,2*N + idopplerP - 20 ):1:min((4*N-1),2*N + idopplerP + 20 ) 
               
            p = ind1 - 2*M ;
            q = ind2 - 2*N ;
            
                % wnu(x,.,.) returns the value of the RRC pulse at x
                f12 = @(tau_prime,nu_prime) (wnu(tau_prime, beta_tau, constbetatau, abstoltau)  .* wnu(tau_prime + delayP(ijk) - p, beta_tau, constbetatau, abstoltau)  .*  wnu(nu_prime, beta_nu, constbetanu, abstolnu)  .*  wnu(nu_prime + dopplerP(ijk) - q, beta_nu, constbetanu, abstolnu)  .* exp(-1i*2*pi*(dopplerP(ijk)*tau_prime/BT))    .* exp(-1i*2*pi*tau_prime.*nu_prime/BT)   .* exp(1i*2*pi*nu_prime*p/BT) );
                hdd(ind1, ind2) = hdd(ind1, ind2) + hgainP(ijk)*exp(-1i*2*pi*delayP(ijk)*dopplerP(ijk)/(BT))* exp(1i*2*pi*p*dopplerP(ijk)/BT) * integral2(f12, -10 + max(0, p - delayP(ijk)), 10 + min(0, p - delayP(ijk) ) ,-10+ max(0, q - dopplerP(ijk)), 10 + min(0, q - dopplerP(ijk) ));
                
          
           end
   
        end
    end
    
    
    
    
    % H_otfs is the effective MN X MN DD domain channel matrix, see the
    % mattrix H_{dd} in equations (2.95) and (2.96) of the book
    % equation (2.96) describes how the elements of H_{d} matrix depend on
    % the effective DD domain channel filter h_{eff}[k,l]    
    H_otfs = zeros(BT, BT);
    for kprime=0:1:(M-1)
        for lprime=0:1:(N-1)
            for k=0:1:(M-1)
                for l=0:1:(N-1)
                    for n=-1:1:1
                        for m=-1:1:1
                    H_otfs(lprime*M + kprime + 1, l*M + k + 1) = H_otfs(lprime*M + kprime + 1, l*M + k + 1) + hdd(2*M + kprime - k - n*M, 2*N+lprime - l - m*N) * exp(1i*2*pi*(lprime - l - m*N)*(k + n*M)/BT) * exp(1i*2*pi*n*l/N);
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
                             H_otfs_est(lprime*M + kprime + 1, l*M + k + 1) = hdd_est(k_index,l_index) *  exp(1i*2*pi*(lprime - l - m*N)*(k + n*M)/BT) * exp(1i*2*pi*n*l/N) ;
                          end
                      end
                  end
               end
            end
        end
    end




    % compute MMSE matrix for equalization at receiver
    MMSE_otfs = inv(H_otfs_est'*H_otfs_est + eye(BT,BT)*(1/SNR_lin)) * H_otfs_est' ;
    
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
    disp('Zak-OTFS BER (RRC)');
    disp(berOTFS/(itr*BT*2));
    

end


berOTFS = berOTFS/(NUM_ITR*BT*2);




   


