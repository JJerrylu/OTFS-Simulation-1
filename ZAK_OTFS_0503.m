clear all;
%----- Program INPUTs---------------------
NUM_ITR = 1; % No. of iterations for BER averaging
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
berOTFS= 0;

for itr=1:1:NUM_ITR  

    for ijk=1:1:numPaths
    dopplerP(ijk) = numax *cos(2*pi*rand);
    end
        
    for pindex=1:1:numPaths
        hgainP(pindex) = sqrt(linpowProfile(pindex))*exp(1i*2*pi*rand);
    end
    
    btau = B*delayP
    tdopplerP = T*dopplerP  
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
    %{
     if itr==1
      figure; bar3(hdd);
      title('H_dd channel matrix 4M*4N'); xlabel('Delay \tau'); ylabel('dopler \nu');
    end
    %}
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
    %{
     if itr==1
      figure; bar3(H_otfs);
      title('H_OTFS channel matrix MN*MN'); xlabel('Delay \tau'); ylabel('dopler \nu');
     end
    %}
    % compute MMSE matrix for equalization at receiver
    MMSE_otfs = inv(H_otfs'*H_otfs + eye(BT,BT)*(1/SNR_lin)) * H_otfs' ;
  
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
    
  
    disp('------------------------');
    disp(itr);
    disp('Zak-OTFS');
    disp(berOTFS/(itr*BT * 2));

end % for itr=1:1:NUM_ITR
    
% Computation of average BER for FDM and Zak-OTFS (perfect CSI)
berOTFS = berOTFS/(NUM_ITR*BT*2);


