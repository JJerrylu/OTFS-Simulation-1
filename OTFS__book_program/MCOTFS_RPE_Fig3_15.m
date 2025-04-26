% =========================================================================
% Book title: OTFS Modulation: Theory and Applications
% Author: Saif Khan Mohammed, Ronny Hadani and Ananthanarayanan
% Chockalingam
% Editor: Wiley-IEEE Press
% Date: Oct. 2024
% Description: This Matlab program generates the heatmap for the Relative
% Prediction Error (RPE) for Multicarrier OTFS (MC-OTFS)as plotted in Fig. 3.15
% of the book for three
% different Zak-OTFS modulations corrsponding to \nu_p = 240 KHz, 30 KHz, 1.25 KHz
% for a two-path channel. For generating the three heatmaps in Fig. 2.27,
% this program needs to be run separately for each of \nu_p = 240 KHz, 30
% KHz, 1.25 KHz.
% RPE is calculated as per equation (2.89) in the book.
% For a given input (i.e., channel delay and Doppler spread,
% bandwidth and time duration of Zak-OTFS pulsone, and \nu_p) the
% program outputs a 2-D heatmap plot of the RPE
% =========================================================================

clear all;

%----- Program INPUTs---------------------
delayspread = 5e-6; % Channel delay spread in seconds
Dopplerspread = 1630; % Channel Doppler spread in Hz
B = 0.96e6; % Pulsone bandwidth in Hz
T = 1.6e-3; % Pulsone time duration in seconds
nup = 240e3; % Zak-OTFS modulation parameter \nu_p in Hz
taup = 1/nup;
%----------------------------------

%------- Channel parameters--------
numPaths = 2; % Number of channel paths is 2
delayP = [ 0 , delayspread ] ;  % The first path induces zero delay shift
dopplerP =  0.5*Dopplerspread*[1 , -1]; % First path induces a Doppler shift of 0.5*Dopplerspread
                                        % second path induces a Doppler shift
                                        % of -0.5*Dopplerspread
hgainP =  1/sqrt(2) * [ 1 , 1]; % Channel gain for each path is 1/sqrt(2)



%------- Zak-OTFS modulation variables
M = round(B*taup); % M = B \tau_p
N = round(T*nup);  % N = T \nu_p
BT = round(B * T) ;
btau = B*delayP; %Normalized path delays
tdopplerP = T*dopplerP; % Normalized path Doppler shift


L = numPaths;
taunorm = delayP/taup;
nunorm = dopplerP/nup;
    
    
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
    

    kp = N/2;
    lp = M/2;
    hdd_est = zeros(N,M);

    for kprime=0:1:(N-1)
        for lprime=0:1:(M-1)
            hdd_est(kprime +1 , lprime+1) = H_otfs(lprime + kprime*M + 1, lp + kp*M + 1);
        end
    end

    H_otfs_est = zeros(M*N, M*N);
    est_error = zeros(M,N);
    for kprime=0:1:(N-1)
        for lprime=0:1:(M-1)
            for k=0:1:(N-1)
                for l=0:1:(M-1)
                    kpt =  mod(kp + kprime -k,N) + 1;
                    lpt =  mod(lp + lprime -l,M) + 1;
                    if ((kpt > 0) && (kpt <= N) && (lpt > 0) && (lpt <= M))
                       H_otfs_est(lprime + kprime*M + 1, l+k*M + 1) = hdd_est(kpt, lpt );
                    end
                               est_error(lprime+1,kprime+1) = est_error(lprime+1, kprime+1) + abs(H_otfs_est(kprime*M + lprime + 1, k*M + l + 1) -  H_otfs(kprime*M + lprime + 1, k*M + l + 1))^2 ;
                end
            end
            
 
            
        end
    end


    
    


      figure(24);
      
      hmtap_hndl = heatmap(10*log10(est_error'/(norm(hdd_est,'fro')^2)),'CellLabelColor','none');
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
      
  
        
    






  


