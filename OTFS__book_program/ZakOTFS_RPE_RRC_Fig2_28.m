clear all;
% =========================================================================
% Book title: OTFS Modulation: Theory and Applications
% Author: Saif Khan Mohammed, Ronny Hadani and Ananthanarayanan
% Chockalingam
% Editor: Wiley-IEEE Press
% Date: Oct. 2024
% Description: This Matlab program generates the heatmap for the Relative
% Prediction Error (RPE) as plotted in Fig 2.27 of the book for three
% different Zak-OTFS modulations corrsponding to \nu_p = 240 KHz, 30 KHz, 1.25 KHz
% for a two-path channel. For generating the three heatmaps in Fig. 2.28,
% this program needs to be run separately for each of \nu_p = 240 KHz, 30
% KHz, 1.25 KHz.
% RPE is calculated as per equation (2.89) in the book. Pulse shaping
% filter is Root Raised Cosine RRC (see equation (2.81) in book).
% For a given input (i.e., channel delay and Doppler spread,
% bandwidth and time duration of Zak-OTFS pulsone, and \nu_p) the
% program outputs a 2-D heatmap plot of the RPE
% =========================================================================


%----- Program INPUTs---------------------
delayspread = 5e-6; % Channel delay spread in seconds
Dopplerspread = 1630; % Channel Doppler spread in Hz
B = 0.96e6; % Pulsone bandwidth in Hz
T = 1.6e-3; % Pulsone time duration in seconds
nup = 30e3; % Zak-OTFS modulation parameter \nu_p (in Hz)
taup = 1/nup;
beta_nu = 0.2; % RRC roll-off factor \beta_nu in book for pulse shape along Doppler domain
beta_tau = 0.1; % RRC roll-off factor \beta_nu in book for pulse shape along delay domain
%----------------------------------

%------- Channel parameters--------
numPaths = 2; % Number of channel paths is 2
delayP = [ 0 , delayspread ]*B ;  % Normalized path delays
dopplerP =  0.5*Dopplerspread*T*[1 , -1]; % Normalized path Doppler shifts
hgainP =  1/sqrt(2) * [ 1 , 1]; % Channel gain for each path is 1/sqrt(2)



%------- Zak-OTFS modulation variables
M = round(B*taup); % M = B \tau_p
N = round(T*nup);  % N = T \nu_p
BT = round(B * T) ;


constbetanu = (beta_nu*sqrt(1/2))* ( (1 + (2/pi))*sin(pi/(4*beta_nu))   +  (1 - (2/pi))*cos(pi/(4*beta_nu))) ;
abstolnu = 0.001;

constbetatau = (beta_tau*sqrt(1/2))* ( (1 + (2/pi))*sin(pi/(4*beta_tau))   +  (1 - (2/pi))*cos(pi/(4*beta_tau))) ;
abstoltau = 0.001;


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
                    % The following equation is based on equation (2.95)
                    % in the book
                    H_otfs(lprime*M + kprime + 1, l*M + k + 1) = H_otfs(lprime*M + kprime + 1, l*M + k + 1) + hdd(2*M + kprime - k - n*M, 2*N+lprime - l - m*N) * exp(1i*2*pi*(lprime - l - m*N)*(k + n*M)/BT) * exp(1i*2*pi*n*l/N);
                        end
                    end
                end
            end
        end
    end



    % pilot location: DD location of pulse transmitted for estimating the taps of h_{eff}[k,l]
    kp = round(M/2);
    lp = round(N/2);

    % Variable hdd_est(.,.) models  \widehat{h}_{eff}[k,l] in the book
    % which is the estimate of the effective DD channel filter h_{eff}
    % based on the received DD channel response to the pilot DD pulse
    % transmitted at (M/2, N/2)

    for kprime=0:1:(M-1)
        for lprime=0:1:(N-1)
           % The following equation implements equation (2.86) in the book
           hdd_est(kp + kprime - kp + 1, lp + lprime -lp + 1) =  H_otfs(lprime*M + kprime + 1, lp*M + kp + 1) * exp(-1i*2*pi*(lprime - lp)*kp/BT);
        end
    end


    % H_otfs_est is the estimated effective DD domain channel matrix based
    % on the estimated channel filter hdd_est (constructed based on
    % equation (2.95) in book)
    H_otfs_est = zeros(BT, round(BT));
    est_error = zeros(M,N);
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
                  % est_error(k^(g)+1, l^(g)+1) is the un-normalized numerator 
                  % in the R.H.S. of equation (2.89) in the book
                  est_error(k+1,l+1) = est_error(k+1, l+1) + abs(H_otfs_est(lprime*M + kprime + 1, l*M + k + 1) -  H_otfs(lprime*M + kprime + 1, l*M + k + 1))^2 ;                  
                  
               end
            end
        end
    end
    


      %Heatmap plot for RPE 


      figure(23);
      % Normalization and plotting the heatmap of RPE (see equation (2.89)
      % in the book)
      hmtap_hndl = heatmap(10*log10(est_error'/(norm(hdd,'fro')^2)),'CellLabelColor','none');
      xv = string((1:M)); %Delay
      yv = string((1:N)); %Doppler
      xv(2:end-1) = "" ;
      yv(2:end-1) = "" ; 
      hmtap_hndl.XDisplayLabels = xv;
      hmtap_hndl.YDisplayLabels = yv;
      hmtap_hndl.XLabel = 'Delay Domain' ;
      hmtap_hndl.YLabel = 'Doppler Domain' ;
      hmtap_hndl.Title = 'Relative Prediction Error (dB)';
      hmtap_hndl.FontSize = 40;
      hmtap_hndl.Colormap = parula;
      hmtap_hndl.ColorLimits = [-80 , 3];
      
  



      

