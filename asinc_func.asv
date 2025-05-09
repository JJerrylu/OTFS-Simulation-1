clc, clear all
T = 1;
delta_f = 1;
M = 14;
N = 12;
nu = 0 :1/(5*N*T) : 1/T-1/(5*N*T);
tau = 0 :1/(10*M*delta_f) : 1/delta_f;

figure;
for n = 1:N
%f_nu1 = sind(rad2deg(pi*nu+ *N*T))./sind(rad2deg(pi*nu*T));
f_nu(n,:) = (1/sqrt(N))*sind(rad2deg(pi*(nu+ (n-1)/(N*T))*N*T))./sind(rad2deg(pi*(nu + (n-1)/(N*T))*T));
hold on
grid on
plot(nu, f_nu(n,:).^2 );
end
hold off
f_nu_sum = sum((f_nu).^2, 1);
figure;
plot(nu, f_nu_sum);
