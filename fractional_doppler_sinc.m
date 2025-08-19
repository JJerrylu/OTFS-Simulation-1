clear; clc; close all;

%% Parameters
M  = 16;              % # subcarriers (delay bins)
N  = 16;              % # time slots (Doppler bins)
Delta_f = 15e3;          % subcarrier spacing (normalized)
T  = 1/Delta_f;       % symbol duration
p  = 0:M-1;
p_ = p;

% kappa values for test
kappa_vals = [0, 0.3];   % fractional Doppler in doppler bin

%% Function: W_{p,p'}(\nu)
W_fun = @(nu) ...
     exp(1j*pi*(nu + (p_'-p)*Delta_f)*T )...
      .* sinc((nu + (p_'-p)*Delta_f)*T);

%% Function: Dirichlet_N
Dirichlet_N = @(x) sin(pi*N*x) ./ (N .* sin(pi*x));

%% Plot |W| for two kappa
figure;
for i = 1:length(kappa_vals)
    nu = kappa_vals(i)./(N*T);  % Doppler in Hz
    Wmat = W_fun(nu);
    subplot(2,2,i);
    imagesc(abs(Wmat));
    colorbar;
    title(sprintf('|W_{p,p''}|, kappa = %.1f', kappa_vals(i)));
    xlabel('p'''); ylabel('p');
end

%% Plot |Dirichlet_N| for two kappa
for i = 1:length(kappa_vals)
    kappa = kappa_vals(i);
    k = 0:N-1;  % Doppler bins
    D = Dirichlet_N(k - kappa);
    % Make it M-by-N just for 2D heatmap visualization
    Dmat = repmat(abs(D), M, 1);
    subplot(2,2,i+2);
    imagesc(Dmat);
    colorbar;
    title(sprintf('|Dirichlet_N(k - kappa)|, kappa = %.1f', kappa_vals(i)));
    xlabel('Doppler bin k'); ylabel('delay bin m');
end
colormap jet;
