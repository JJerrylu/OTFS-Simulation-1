clc, clear all
M = 64;
N = 64;
MN = M*N;
delta_f = 15e3;            % 15kHz
T = 1/delta_f;             % Block Duration
fc = 4e9;
delay_tap = [0,1,2,3,4,5,8,13,19];
tap = length(delay_tap);
%% Delay and channel coefficient
delay = delay_tap/(M*delta_f);
delay_pow_prof = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];
lin_pow_prof = 10.^(delay_pow_prof/10);
lin_pow_prof = lin_pow_prof/sum(lin_pow_prof);
chan_coef = sqrt(lin_pow_prof).*(sqrt(1/2) * (randn(1,tap)+1i*randn(1,tap)));
%% Doppler shift
max_speed = 500;           % 500 km/h
v = max_speed*(1000/3600);
max_doppler = (v*fc)/(299792458);
doppler = max_doppler*cos(2*pi*rand(1,tap));
doppler_inx = doppler*(N*T); % Doppler taps using Jake's spectrum
%doppler_inx = [1, 2, -4, -3, 5, 6, -10, 13, -14];
%doppler = doppler_inx/(N*T);
int_doppler_inx = round(doppler_inx);                  % Integer Doppler
frac_doppler_inx = doppler_inx-int_doppler_inx;        % Fractional Doppler

X  = zeros(M, N);
H  = zeros(M, N);

X(M/2, N/2) = 1;
%% rectangular pulse shaping DD domain pilot response
Ni = 16;
for l = 0:M-1
    for k = 0:N-1
        for i = 1:tap
            %% if integer doppler 
            if frac_doppler_inx(i) == 0
                      H(l+1, k+1) = H(l+1,k+1)+ chan_coef(i)*exp(-1i*2*pi*(l-delay_tap(i))*doppler(i))...
                    * X(mod(l-delay_tap(i), M)+1, mod(k-int_doppler_inx(i), N)+1);
                      continue;
            end
            %% fractional doppler
            for q =  -Ni:Ni
                beta = (exp(-1i*2*pi*(-q-frac_doppler_inx(i)))-1)/(exp(-1i*2*pi/N*(-q-frac_doppler_inx(i)))-1);
                if l>=delay_tap(i) && l<M
                    H(l+1, k+1) = H(l+1,k+1)+ chan_coef(i)*exp(-1i*2*pi*(l-delay_tap(i))*doppler(i))...
                    * 1/N*(beta)...
                    * X(mod(l-delay_tap(i), M)+1, mod(k-int_doppler_inx(i)+q, N)+1);
                else
                    H(l+1, k+1) = H(l+1,k+1)+ chan_coef(i)*exp(-1i*2*pi*(l-delay_tap(i))*doppler(i))...
                    * 1/N*(beta-1)*exp(-1i*2*pi*mod((k+1-int_doppler_inx(i)+q),N))...
                    * X(mod(l-delay_tap(i), M)+1, mod(k-int_doppler_inx(i)+q, N)+1);
                end
            end
        
        end
    end
end
figure; bar3(abs(X));
figure; bar3(abs(H)); title('Delay Doppler domain channel matrix H\_dd (abs)'); xlabel('doppler'); ylabel('delay');
