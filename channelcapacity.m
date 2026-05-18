
%% Channel Capacity and Outage Probability Analysis
% Models: AWGN, Rayleigh, Rician, and Eta-Mu
clear; clc; close all;

%% Parameters
SNR_dB = 0:2:30;             % SNR range in dB
SNR_lin = 10.^(SNR_dB/10);   % Linear SNR
N = 10^5;                    % Number of Monte Carlo iterations
R_target = 2;                % Target rate for outage (bits/s/Hz)

% Channel Specific Parameters
K_rician = 5;                % Rician K-factor (ratio of LOS to NLOS)
eta = 0.5;                   % Eta-mu: Power ratio of I and Q components
mu = 1.5;                    % Eta-mu: Number of multipath clusters

%% Pre-allocate Arrays
cap_awgn = zeros(size(SNR_dB));
cap_rayleigh = zeros(size(SNR_dB));
cap_rician = zeros(size(SNR_dB));
cap_etamu = zeros(size(SNR_dB));

out_awgn = zeros(size(SNR_dB));
out_rayleigh = zeros(size(SNR_dB));
out_rician = zeros(size(SNR_dB));
out_etamu = zeros(size(SNR_dB));

%% Monte Carlo Simulation
fprintf('Starting Monte Carlo Simulation...\n');

for i = 1:length(SNR_lin)
    snr = SNR_lin(i);
    
    % 1. AWGN Channel (H = 1 constant)
    h_awgn = ones(1, N);
    C_awgn = log2(1 + snr * abs(h_awgn).^2);
    cap_awgn(i) = mean(C_awgn);
    out_awgn(i) = sum(C_awgn < R_target) / N;

    % 2. Rayleigh Fading (Complex Gaussian)
    % E[|h|^2] = 1
    h_rayleigh = (randn(1, N) + 1i*randn(1, N)) / sqrt(2);
    C_rayleigh = log2(1 + snr * abs(h_rayleigh).^2);
    cap_rayleigh(i) = mean(C_rayleigh);
    out_rayleigh(i) = sum(C_rayleigh < R_target) / N;

    % 3. Rician Fading
    % s^2 is LOS power, sigma^2 is scattered power. K = s^2 / (2*sigma^2)
    s = sqrt(K_rician / (K_rician + 1));
    sigma = sqrt(1 / (2 * (K_rician + 1)));
    h_rician = (s + sigma*randn(1,N)) + 1i*(sigma*randn(1,N));
    C_rician = log2(1 + snr * abs(h_rician).^2);
    cap_rician(i) = mean(C_rician);
    out_rician(i) = sum(C_rician < R_target) / N;

    % 4. Eta-Mu Fading (Format 1)
    % Generate using the sum of squared Gaussian components
    % Power is distributed between 2*mu channels with variance ratio eta
    sigma_i = sqrt(1 / (2 * mu * (1 + eta))); % Variance of In-phase
    sigma_q = sqrt(eta / (2 * mu * (1 + eta))); % Variance of Quadrature
    
    H_etamu_sq = zeros(1, N);
    for m = 1:ceil(mu) % Simplified realization of eta-mu
        X = sigma_i * randn(1, N);
        Y = sigma_q * randn(1, N);
        H_etamu_sq = H_etamu_sq + X.^2 + Y.^2;
    end
    % Normalize power to 1 (accounting for integer mu approximation)
    H_etamu_sq = H_etamu_sq / mean(H_etamu_sq); 
    
    C_etamu = log2(1 + snr * H_etamu_sq);
    cap_etamu(i) = mean(C_etamu);
    out_etamu(i) = sum(C_etamu < R_target) / N;
end

%% Plotting Results
figure('Position', [100, 100, 1000, 450]);

% Subplot 1: Ergodic Capacity
subplot(1, 2, 1);
plot(SNR_dB, cap_awgn, 'k--+', 'LineWidth', 1.5); hold on;
plot(SNR_dB, cap_rayleigh, 'r-o', 'LineWidth', 1.5);
plot(SNR_dB, cap_rician, 'b-s', 'LineWidth', 1.5);
plot(SNR_dB, cap_etamu, 'm-d', 'LineWidth', 1.5);
grid on;
xlabel('Average SNR (dB)');
ylabel('Capacity (bits/s/Hz)');
title('Ergodic Channel Capacity');
legend('AWGN (Shannon Limit)', 'Rayleigh', 'Rician (K=5)', '\eta-\mu (\eta=0.5, \mu=1.5)', 'Location', 'NorthWest');

% Subplot 2: Outage Probability
subplot(1, 2, 2);
semilogy(SNR_dB, out_awgn, 'k--+', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB, out_rayleigh, 'r-o', 'LineWidth', 1.5);
semilogy(SNR_dB, out_rician, 'b-s', 'LineWidth', 1.5);
semilogy(SNR_dB, out_etamu, 'm-d', 'LineWidth', 1.5);
grid on;
ylim([1e-4 1]);
xlabel('Average SNR (dB)');
ylabel('Outage Probability P(C < R)');
title(['Outage Probability (Target Rate R = ', num23str(R_target), ')']);
legend('AWGN', 'Rayleigh', 'Rician', '\eta-\mu', 'Location', 'SouthWest');

fprintf('Simulation Complete.\n');