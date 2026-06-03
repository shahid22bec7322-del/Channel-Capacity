%% ========================================================================
%  Channel Capacity Analysis under Fading Channels
%  Compares: AWGN, Rayleigh, Rician, and Eta-Mu Fading Channels
%  Methods: Monte Carlo Simulation + Analytical Bounds
%  Author: Channel Analysis Suite
%% ========================================================================

clc; clear; close all;

%% ========================== PARAMETERS ==================================
SNR_dB       = -10:2:30;          % SNR range in dB
SNR_lin      = 10.^(SNR_dB/10);   % Linear SNR
N_mc         = 1e5;               % Monte Carlo iterations
outage_thresh_dB = 0;             % Outage threshold SNR (dB) => C_th = log2(1+1)=1 bit/s/Hz
C_threshold  = log2(1 + 10^(outage_thresh_dB/10)); % bits/s/Hz

% Rician K-factor (dB -> linear)
K_rice       = 3;                  % K = 3 (moderate LOS component)

% Eta-Mu parameters
eta_param    = 0.5;                % Eta: ratio of in-phase to quadrature variances
mu_param     = 1.5;                % Mu: number of multipath clusters

fprintf('=== Channel Capacity Analysis ===\n');
fprintf('SNR Range   : %d to %d dB\n', SNR_dB(1), SNR_dB(end));
fprintf('Monte Carlo : %d samples\n', N_mc);
fprintf('Outage C_th : %.4f bits/s/Hz\n\n', C_threshold);

%% =================== PREALLOCATE RESULT ARRAYS ==========================
n_snr = length(SNR_dB);

C_awgn      = zeros(1, n_snr);
C_rayleigh  = zeros(1, n_snr);
C_rician    = zeros(1, n_snr);
C_etamu     = zeros(1, n_snr);

P_out_awgn     = zeros(1, n_snr);
P_out_rayleigh = zeros(1, n_snr);
P_out_rician   = zeros(1, n_snr);
P_out_etamu    = zeros(1, n_snr);

%% =================== MONTE CARLO SIMULATION =============================
fprintf('Running Monte Carlo simulations...\n');

for i = 1:n_snr
    snr = SNR_lin(i);

    %% --- 1. AWGN Channel ---
    % Deterministic: no fading, capacity = log2(1 + SNR)
    C_awgn(i)      = log2(1 + snr);
    P_out_awgn(i)  = (C_awgn(i) < C_threshold); % binary: 0 or 1

    %% --- 2. Rayleigh Fading Channel ---
    % |h|^2 ~ Exp(1), instantaneous SNR gamma = snr * |h|^2
    h_rayleigh     = (1/sqrt(2)) * (randn(1,N_mc) + 1j*randn(1,N_mc));
    gamma_ray      = snr * abs(h_rayleigh).^2;
    C_inst_ray     = log2(1 + gamma_ray);
    C_rayleigh(i)  = mean(C_inst_ray);
    P_out_rayleigh(i) = mean(C_inst_ray < C_threshold);

    %% --- 3. Rician Fading Channel ---
    % LOS component s = sqrt(K/(K+1)), scatter sigma = sqrt(1/(2*(K+1)))
    s_rice  = sqrt(K_rice / (K_rice + 1));
    sig_rice = sqrt(1 / (2*(K_rice + 1)));
    h_rician = (s_rice + sig_rice*randn(1,N_mc)) + 1j*(sig_rice*randn(1,N_mc));
    gamma_ric      = snr * abs(h_rician).^2;
    C_inst_ric     = log2(1 + gamma_ric);
    C_rician(i)    = mean(C_inst_ric);
    P_out_rician(i) = mean(C_inst_ric < C_threshold);

    %% --- 4. Eta-Mu Fading Channel ---
    % Eta-Mu: generalized fading model. Samples generated via:
    %   X = sqrt(mu/(2*h)) * (chi_2mu_I + eta*chi_2mu_Q) 
    % where h = (2+eta+1/eta)*mu/4
    % Envelope squared follows Eta-Mu PDF.
    %
    % Efficient generation: sum of scaled chi-squared RVs
    %   In-phase  power: P_I ~ Gamma(mu, 1/(2*mu*h_I)) where h_I = mu/(1+eta)
    %   Quad power: P_Q ~ Gamma(mu, 1/(2*mu*h_Q)) where h_Q = mu*eta/(1+eta)
    %
    % Using format 1 of Eta-Mu:
    %   h = (2+eta+1/eta)*mu/4
    %   Omega = E[r^2] = 1 (normalized)
    %   P_I ~ Gamma(mu, (1+eta)/(2*mu))
    %   P_Q ~ Gamma(mu, (1+1/eta)/(2*mu))

    shape_I  = mu_param;
    scale_I  = (1 + eta_param) / (2 * mu_param);
    shape_Q  = mu_param;
    scale_Q  = (1 + 1/eta_param) / (2 * mu_param);

    P_I = gamrnd(shape_I, scale_I, 1, N_mc);
    P_Q = gamrnd(shape_Q, scale_Q, 1, N_mc);

    % Envelope squared = P_I + P_Q (normalized power)
    r2_etamu       = P_I + P_Q;
    gamma_etamu    = snr * r2_etamu;
    C_inst_etamu   = log2(1 + gamma_etamu);
    C_etamu(i)     = mean(C_inst_etamu);
    P_out_etamu(i) = mean(C_inst_etamu < C_threshold);

    if mod(i, 5) == 0
        fprintf('  Processed SNR = %+d dB (%d/%d)\n', SNR_dB(i), i, n_snr);
    end
end

fprintf('\nSimulation complete.\n\n');

%% =================== DIVERSITY & CAPACITY LOSS ANALYSIS =================
fprintf('=== Diversity Gain & Capacity Loss (at SNR = 20 dB) ===\n');
idx_ref = find(SNR_dB == 20);
if isempty(idx_ref), idx_ref = n_snr; end

cap_loss_ray = C_awgn(idx_ref) - C_rayleigh(idx_ref);
cap_loss_ric = C_awgn(idx_ref) - C_rician(idx_ref);
cap_loss_eta = C_awgn(idx_ref) - C_etamu(idx_ref);

fprintf('  AWGN Capacity      : %.4f bits/s/Hz\n', C_awgn(idx_ref));
fprintf('  Rayleigh Capacity  : %.4f bits/s/Hz (Loss: %.4f)\n', C_rayleigh(idx_ref), cap_loss_ray);
fprintf('  Rician Capacity    : %.4f bits/s/Hz (Loss: %.4f)\n', C_rician(idx_ref), cap_loss_ric);
fprintf('  Eta-Mu Capacity    : %.4f bits/s/Hz (Loss: %.4f)\n\n', C_etamu(idx_ref), cap_loss_eta);

% Outage probability at 10 dB
idx_10 = find(SNR_dB == 10);
if isempty(idx_10), idx_10 = round(n_snr/2); end
fprintf('=== Outage Probability at SNR = 10 dB ===\n');
fprintf('  Rayleigh  : %.6f (%.2f%%)\n', P_out_rayleigh(idx_10), P_out_rayleigh(idx_10)*100);
fprintf('  Rician    : %.6f (%.2f%%)\n', P_out_rician(idx_10), P_out_rician(idx_10)*100);
fprintf('  Eta-Mu    : %.6f (%.2f%%)\n\n', P_out_etamu(idx_10), P_out_etamu(idx_10)*100);

%% ========================== ANALYTICAL CHECKS ===========================
% Rayleigh analytical outage: P_out = 1 - exp(-gamma_th/SNR)
gamma_th_outage = 2^C_threshold - 1;
P_out_ray_analytic = 1 - exp(-gamma_th_outage ./ SNR_lin);

%% ========================== PLOTTING ====================================
colors = struct(...
    'awgn',    [0.12 0.47 0.71], ...   % Blue
    'rayleigh',[0.89 0.10 0.11], ...   % Red
    'rician',  [0.20 0.63 0.17], ...   % Green
    'etamu',   [1.00 0.50 0.00]);       % Orange

fig_pos = [100 100 1400 900];

%% --- FIGURE 1: Ergodic Channel Capacity ---
figure('Position', fig_pos, 'Name', 'Ergodic Channel Capacity');
plot(SNR_dB, C_awgn,     '-o', 'Color', colors.awgn,    'LineWidth',2.5, 'MarkerSize',5, 'DisplayName','AWGN (Deterministic)');
hold on;
plot(SNR_dB, C_rayleigh, '-s', 'Color', colors.rayleigh, 'LineWidth',2.5, 'MarkerSize',5, 'DisplayName','Rayleigh (K=0)');
plot(SNR_dB, C_rician,   '-^', 'Color', colors.rician,  'LineWidth',2.5, 'MarkerSize',5, 'DisplayName',sprintf('Rician (K=%d)', K_rice));
plot(SNR_dB, C_etamu,    '-d', 'Color', colors.etamu,   'LineWidth',2.5, 'MarkerSize',5, 'DisplayName',sprintf('\\eta-\\mu (\\eta=%.1f, \\mu=%.1f)', eta_param, mu_param));
yline(C_threshold, '--k', 'LineWidth', 1.5, 'DisplayName', sprintf('Outage Threshold (%.2f b/s/Hz)', C_threshold));
grid on; box on;
xlabel('Average SNR (dB)', 'FontSize', 14);
ylabel('Ergodic Capacity (bits/s/Hz)', 'FontSize', 14);
title('Ergodic Channel Capacity Comparison', 'FontSize', 15, 'FontWeight','bold');
legend('Location','northwest', 'FontSize', 11);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
xlim([SNR_dB(1) SNR_dB(end)]);
ylim([0 max(C_awgn)*1.05]);

%% --- FIGURE 2: Outage Probability ---
figure('Position', fig_pos, 'Name', 'Outage Probability');
semilogy(SNR_dB, max(P_out_rayleigh, 1e-6), '-s', 'Color', colors.rayleigh, 'LineWidth',2.5, 'MarkerSize',5, 'DisplayName','Rayleigh (MC)');
hold on;
semilogy(SNR_dB, max(P_out_ray_analytic, 1e-6), '--', 'Color', colors.rayleigh, 'LineWidth',1.5, 'DisplayName','Rayleigh (Analytical)');
semilogy(SNR_dB, max(P_out_rician, 1e-6),   '-^', 'Color', colors.rician,  'LineWidth',2.5, 'MarkerSize',5, 'DisplayName',sprintf('Rician (K=%d)', K_rice));
semilogy(SNR_dB, max(P_out_etamu, 1e-6),    '-d', 'Color', colors.etamu,   'LineWidth',2.5, 'MarkerSize',5, 'DisplayName',sprintf('\\eta-\\mu (\\eta=%.1f, \\mu=%.1f)', eta_param, mu_param));
grid on; box on;
xlabel('Average SNR (dB)', 'FontSize', 14);
ylabel('Outage Probability P_{out}', 'FontSize', 14);
title('Outage Probability vs. SNR', 'FontSize', 15, 'FontWeight','bold');
legend('Location','southwest', 'FontSize', 11);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
xlim([SNR_dB(1) SNR_dB(end)]);
ylim([1e-5 1]);

%% --- FIGURE 3: Capacity Loss due to Fading ---
figure('Position', fig_pos, 'Name', 'Capacity Loss');
loss_ray = C_awgn - C_rayleigh;
loss_ric = C_awgn - C_rician;
loss_eta = C_awgn - C_etamu;

plot(SNR_dB, loss_ray, '-s', 'Color', colors.rayleigh, 'LineWidth',2.5, 'MarkerSize',5, 'DisplayName','Rayleigh Loss');
hold on;
plot(SNR_dB, loss_ric, '-^', 'Color', colors.rician,  'LineWidth',2.5, 'MarkerSize',5, 'DisplayName',sprintf('Rician Loss (K=%d)', K_rice));
plot(SNR_dB, loss_eta, '-d', 'Color', colors.etamu,   'LineWidth',2.5, 'MarkerSize',5, 'DisplayName',sprintf('\\eta-\\mu Loss'));
yline(0, '-k', 'LineWidth', 1.2);
grid on; box on;
xlabel('Average SNR (dB)', 'FontSize', 14);
ylabel('Capacity Loss vs. AWGN (bits/s/Hz)', 'FontSize', 14);
title('Capacity Loss Due to Fading', 'FontSize', 15, 'FontWeight','bold');
legend('Location','northwest', 'FontSize', 11);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
xlim([SNR_dB(1) SNR_dB(end)]);

%% =================== DIVERSITY GAIN COMPUTATION =========================
% Diversity gain estimated from slope of log(P_out) vs log(SNR) at high SNR
% d = -d(log P_out)/d(log SNR) 
log_snr = log10(SNR_lin(end-5:end));

div_ray = -mean(diff(log10(max(P_out_rayleigh(end-5:end),1e-10))) ./ diff(log_snr));
div_ric = -mean(diff(log10(max(P_out_rician(end-5:end),1e-10))) ./ diff(log_snr));
div_eta = -mean(diff(log10(max(P_out_etamu(end-5:end),1e-10))) ./ diff(log_snr));

fprintf('=== Estimated Diversity Gain (high-SNR slope) ===\n');
fprintf('  Rayleigh : %.3f  (theoretical: 1)\n', div_ray);
fprintf('  Rician   : %.3f  (approaches 1 for K<inf)\n', div_ric);
fprintf('  Eta-Mu   : %.3f  (theoretical: mu = %.1f)\n', div_eta, mu_param);

%% =================== SAVE RESULTS =======================================
results.SNR_dB          = SNR_dB;
results.C_awgn          = C_awgn;
results.C_rayleigh      = C_rayleigh;
results.C_rician        = C_rician;
results.C_etamu         = C_etamu;
results.P_out_rayleigh  = P_out_rayleigh;
results.P_out_rician    = P_out_rician;
results.P_out_etamu     = P_out_etamu;
results.parameters.N_mc       = N_mc;
results.parameters.K_rice     = K_rice;
results.parameters.eta        = eta_param;
results.parameters.mu         = mu_param;
results.parameters.C_threshold = C_threshold;

save('channel_capacity_results.mat', 'results');
fprintf('\nResults saved to: channel_capacity_results.mat\n');
fprintf('All plots rendered successfully.\n');
