clc, clear all
fc =4e9;
delta_f = 15e3;
T = 1/delta_f;
M = 14;
N = 64;
nu = 0 :1/(N*T) : (N-1)/(N*T);
tau = 0 :1/(10*M*delta_f) : 1/delta_f;
tap = 3;
        % Fractional Doppler

%{
for n = 1:N
%f_nu1 = sind(rad2deg(pi*nu*N*T))./sind(rad2deg(pi*nu*T));
f_nu(n,:) = (1/sqrt(N))*sind(rad2deg(pi*(nu+ (n-1)/(N*T))*N*T))./sind(rad2deg(pi*(nu + (n-1)/(N*T))*T));
hold on
grid on
plot(nu, f_nu(n,:).^2 );
end
hold off
f_nu_sum = sum((f_nu).^2, 1);
figure;
plot(nu, f_nu_sum);
%}
x = zeros(1, N);
x(N/2) = 1;
doppler = zeros(1,tap);
doppler_inx = zeros(1, tap);
int_doppler_inx = zeros(1,tap);
frac_doppler_inx = zeros(1,tap);
chan_coef = zeros(1,tap);
beta = zeros(3, N);
%% channel parameter
delay_pow_prof = [0, -1.5, -2];
lin_pow_prof = 10.^(delay_pow_prof/10);
lin_pow_prof = lin_pow_prof/sum(lin_pow_prof);
max_speed = 500;           % 500 km/h
v = max_speed*(1000/3600);
max_doppler = (v*fc)/(299792458);

grid on;
subplot(311)
stem(nu, abs(x), "filled");

subplot(312)
hold on;
for i = 1:3
%% channel coef
chan_coef(i) = sqrt(lin_pow_prof(i)).*(sqrt(1/2) * (randn(1)+1i*randn(1)));
%% Doppler shift
doppler(i) = max_doppler*cos(2*pi*rand(1));
doppler_inx(i) = doppler(i)*(N*T); % Doppler taps using Jake's spectrum
%doppler_inx = [1, 2, -4, -3, 5, 6, -10, 13, -14];
%doppler = doppler_inx/(N*T);
int_doppler_inx(i) = round(doppler_inx(i));                  % Integer Doppler
frac_doppler_inx(i) = doppler_inx(i)-int_doppler_inx(i);
for k = 0:N-1
        if  frac_doppler_inx(i) ==0
            beta(i, k+1) = beta(i,k+1) + chan_coef(i)*x(mod(k-int_doppler_inx(i), N)+1);
            continue;
        end
        for q = 0:N-1
            beta(i, k+1) = beta(i,k+1) + chan_coef(i)/N*(exp(-1i*2*pi*(-q-frac_doppler_inx(i)))-1)/(exp(-1i*2*pi/N*(-q-frac_doppler_inx(i)))-1)...
                        *x(mod(k-int_doppler_inx(i)+q, N)+1);
%beta(i+1) = exp(-1i*2*pi*mod((k-int_doppler_inx(i)+q),N))
        end
end
stem(nu, abs(beta(i,:)), "filled");
end
hold off;

total_path_response = sum(beta,1);
total_path_energy = sum(abs(total_path_response).^2);
subplot(313)
stem(nu, abs(total_path_response), "filled");
