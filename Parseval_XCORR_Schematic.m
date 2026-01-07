%% Band-Partitioned Cross-Correlation: Methodology Demonstration
% Zein Sadek
%
% This script demonstrates the mathematical framework for computing
% band-partitioned cross-correlation coefficients using Parseval's theorem
% and the Wiener-Khinchin relationship.
%
% Key relationships:
%   Parseval:        Var(x) = ∫ Sxx(f) df
%   Wiener-Khinchin: Cov(x,y) = ∫ Re{Sxy(f)} df
%   Correlation:     ρ = Cov(x,y) / sqrt(Var(x) * Var(y))

clear; close all; clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE DEMONSTRATIVE SIGNALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signal parameters
fs = 30;            % Sampling frequency [Hz]
T = 120;            % Duration [seconds]
t = (0:1/fs:T-1/fs)';
n = length(t);

% Create two correlated signals with distinct frequency content
% Signal A: Low-frequency drift + wave-frequency oscillation + noise
% Signal B: Correlated wave component (phase-shifted) + different LF + noise

rng(42);  % For reproducibility

% Frequency components
f_LF = 0.15;        % Low-frequency component [Hz]
f_wave = 1.4;       % Wave frequency [Hz] (similar to LM5)
f_HF = 4.5;         % High-frequency component [Hz]

% Signal A: Platform motion of upstream turbine
A_LF = 0.8 * sin(2*pi*f_LF*t + 0.3);                    % Slow drift
A_wave = 2.0 * sin(2*pi*f_wave*t);                       % Wave response
A_HF = 0.3 * sin(2*pi*f_HF*t + 1.2);                    % High-freq vibration
A_noise = 0.4 * randn(n, 1);                             % Broadband noise
A = A_LF + A_wave + A_HF + A_noise;

% Signal B: Platform motion of downstream turbine (correlated via waves)
phase_lag = 0.4;    % Phase lag due to wave propagation [rad]
B_LF = 0.5 * sin(2*pi*f_LF*t + 1.8);                    % Different LF (uncorrelated)
B_wave = 1.8 * sin(2*pi*f_wave*t - phase_lag);          % Correlated wave response
B_HF = 0.25 * sin(2*pi*f_HF*t + 0.5);                   % Partially correlated HF
B_noise = 0.4 * randn(n, 1);                             % Independent noise
B = B_LF + B_wave + B_HF + B_noise;

% Remove mean (required for correlation = covariance relationship)
A = A - mean(A);
B = B - mean(B);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE FREQUENCY BANDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bands.LF = [0, 0.5];           % Low-frequency band
bands.WF = [0.5, 3.0];         % Wave-frequency band
bands.HF = [3.0, fs/2];        % High-frequency band

band_names = fieldnames(bands);
num_bands = length(band_names);
band_colors = {'#3498db', '#e74c3c', '#2ecc71'};
band_labels = {'Low Frequency (LF)', 'Wave Frequency (WF)', 'High Frequency (HF)'};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME-DOMAIN STATISTICS (Ground Truth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_A = var(A);
var_B = var(B);
cov_AB = mean(A .* B);  % Since mean = 0, this equals covariance
rho_AB = cov_AB / sqrt(var_A * var_B);

fprintf('=== TIME-DOMAIN STATISTICS ===\n')
fprintf('Var(A):   %.4f\n', var_A)
fprintf('Var(B):   %.4f\n', var_B)
fprintf('Cov(A,B): %.4f\n', cov_AB)
fprintf('ρ(A,B):   %.4f\n\n', rho_AB)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE CROSS-CORRELATION FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxLag_s = 5;  % seconds
maxLag = round(maxLag_s * fs);
[r_AB, lags] = xcorr(A, B, maxLag, 'coeff');
tau = lags / fs;

% Find peak
[~, idx_peak] = max(abs(r_AB));
tau_peak = tau(idx_peak);
rho_peak = r_AB(idx_peak);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE SPECTRA (Manual FFT - Parseval Consistent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfft = 2^nextpow2(n);

% FFT of signals
X = fft(A, nfft);
Y = fft(B, nfft);

% Two-sided spectral densities
Sxx_2s = (abs(X).^2) / (fs * n);
Syy_2s = (abs(Y).^2) / (fs * n);
Sxy_2s = (X .* conj(Y)) / (fs * n);

% Frequency vector
f_2s = (0:nfft-1)' * fs / nfft;

% Convert to one-sided
n_os = nfft/2 + 1;
f = f_2s(1:n_os);
Sxx = Sxx_2s(1:n_os);
Syy = Syy_2s(1:n_os);
Sxy = Sxy_2s(1:n_os);

% Double non-DC, non-Nyquist for one-sided
Sxx(2:end-1) = 2 * Sxx(2:end-1);
Syy(2:end-1) = 2 * Syy(2:end-1);
Sxy(2:end-1) = 2 * Sxy(2:end-1);

df = fs / nfft;

% Coherence and phase
Cxy = abs(Sxy).^2 ./ (Sxx .* Syy);
Cxy(isnan(Cxy) | isinf(Cxy)) = 0;
phi_xy = angle(Sxy);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERIFY PARSEVAL'S THEOREM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_A_freq = sum(Sxx) * df;
var_B_freq = sum(Syy) * df;
cov_AB_freq = sum(real(Sxy)) * df;
rho_AB_freq = cov_AB_freq / sqrt(var_A_freq * var_B_freq);

fprintf('=== PARSEVAL VERIFICATION ===\n')
fprintf('Var(A):   Time=%.4f, Freq=%.4f, Error=%.2f%%\n', ...
    var_A, var_A_freq, 100*abs(var_A_freq - var_A)/var_A)
fprintf('Var(B):   Time=%.4f, Freq=%.4f, Error=%.2f%%\n', ...
    var_B, var_B_freq, 100*abs(var_B_freq - var_B)/var_B)
fprintf('Cov(A,B): Time=%.4f, Freq=%.4f, Error=%.2f%%\n', ...
    cov_AB, cov_AB_freq, 100*abs(cov_AB_freq - cov_AB)/abs(cov_AB))
fprintf('ρ(A,B):   Time=%.4f, Freq=%.4f\n\n', rho_AB, rho_AB_freq)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE BAND-PARTITIONED QUANTITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b = 1:num_bands
    bname = band_names{b};
    f_lo = bands.(bname)(1);
    f_hi = bands.(bname)(2);
    
    idx = (f >= f_lo) & (f < f_hi);
    
    % Band-integrated quantities
    band_results.(bname).var_A = sum(Sxx(idx)) * df;
    band_results.(bname).var_B = sum(Syy(idx)) * df;
    band_results.(bname).cov_AB = sum(real(Sxy(idx))) * df;
    
    % Band correlation coefficient
    denom = sqrt(band_results.(bname).var_A * band_results.(bname).var_B);
    band_results.(bname).rho = band_results.(bname).cov_AB / denom;
    
    % Fractions of total
    band_results.(bname).frac_var_A = band_results.(bname).var_A / var_A_freq;
    band_results.(bname).frac_var_B = band_results.(bname).var_B / var_B_freq;
    band_results.(bname).frac_cov = band_results.(bname).cov_AB / cov_AB_freq;
end

% Print band results
fprintf('=== BAND-PARTITIONED RESULTS ===\n')
fprintf('%-6s %10s %10s %10s %10s %10s %10s %10s\n', ...
    'Band', 'Var(A)', 'Var(B)', 'Cov(A,B)', 'ρ_band', '%Var(A)', '%Var(B)', '%Cov')
fprintf('%s\n', repmat('-', 1, 86))
for b = 1:num_bands
    bname = band_names{b};
    fprintf('%-6s %10.4f %10.4f %10.4f %+10.4f %9.1f%% %9.1f%% %9.1f%%\n', ...
        bname, ...
        band_results.(bname).var_A, ...
        band_results.(bname).var_B, ...
        band_results.(bname).cov_AB, ...
        band_results.(bname).rho, ...
        100*band_results.(bname).frac_var_A, ...
        100*band_results.(bname).frac_var_B, ...
        100*band_results.(bname).frac_cov)
end
fprintf('\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1: TIME-DOMAIN SIGNALS AND CROSS-CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'white', 'Position', [50, 50, 1400, 900])
sgtitle('\textbf{Band-Partitioned Cross-Correlation: Time-Domain View}', ...
    'Interpreter', 'latex', 'FontSize', 16)

% Panel 1: Signal A
subplot(3, 2, 1)
plot(t, A, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 0.8)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('$A(t)$', 'Interpreter', 'latex')
title('Signal A (Upstream Turbine)', 'Interpreter', 'latex')
xlim([0, 30])
grid on
box on

% Panel 2: Signal B
subplot(3, 2, 2)
plot(t, B, 'Color', [0.8, 0.3, 0.2], 'LineWidth', 0.8)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('$B(t)$', 'Interpreter', 'latex')
title('Signal B (Downstream Turbine)', 'Interpreter', 'latex')
xlim([0, 30])
grid on
box on

% Panel 3: Overlay of both signals (zoomed)
subplot(3, 2, 3)
plot(t, A, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 1.2, 'DisplayName', 'Signal A')
hold on
plot(t, B, 'Color', [0.8, 0.3, 0.2], 'LineWidth', 1.2, 'DisplayName', 'Signal B')
hold off
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('Signals Overlay (Zoomed)', 'Interpreter', 'latex')
xlim([10, 20])
legend('Interpreter', 'latex', 'Location', 'northeast')
grid on
box on

% Panel 4: Scatter plot A vs B
subplot(3, 2, 4)
scatter(A, B, 10, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.3)
hold on
% Add regression line
p = polyfit(A, B, 1);
A_fit = linspace(min(A), max(A), 100);
B_fit = polyval(p, A_fit);
plot(A_fit, B_fit, 'r-', 'LineWidth', 2)
hold off
xlabel('$A(t)$', 'Interpreter', 'latex')
ylabel('$B(t)$', 'Interpreter', 'latex')
title(sprintf('Scatter Plot: $\\rho = %.3f$', rho_AB), 'Interpreter', 'latex')
axis equal
grid on
box on

% Panel 5: Cross-correlation function
subplot(3, 2, [5, 6])
plot(tau, r_AB, 'k-', 'LineWidth', 1.5)
hold on
plot(tau_peak, rho_peak, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
xline(0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1)
yline(0, '-', 'Color', [0.7, 0.7, 0.7])
hold off
xlabel('Time Lag $\tau$ [s]', 'Interpreter', 'latex')
ylabel('$\rho_{AB}(\tau)$', 'Interpreter', 'latex')
title(sprintf('Cross-Correlation Function: Peak $\\rho = %.3f$ at $\\tau = %.2f$ s', ...
    rho_peak, tau_peak), 'Interpreter', 'latex')
xlim([-maxLag_s, maxLag_s])
ylim([-0.5, 1])
grid on
box on

% Add annotation
text(0.02, 0.95, sprintf('Zero-lag correlation: $\\rho_0 = %.3f$', r_AB(lags == 0)), ...
    'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', 11, ...
    'BackgroundColor', 'white', 'EdgeColor', 'k')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2: FREQUENCY-DOMAIN ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'white', 'Position', [100, 50, 1400, 900])
sgtitle('\textbf{Band-Partitioned Cross-Correlation: Frequency-Domain View}', ...
    'Interpreter', 'latex', 'FontSize', 16)

% Panel 1: PSD of Signal A
subplot(2, 3, 1)
semilogy(f, Sxx, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 1.2)
hold on
for b = 1:num_bands
    idx = (f >= bands.(band_names{b})(1)) & (f < bands.(band_names{b})(2));
    area(f(idx), Sxx(idx), 'FaceColor', band_colors{b}, 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName', band_labels{b})
end
hold off
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('$S_{AA}(f)$ [units$^2$/Hz]', 'Interpreter', 'latex')
title(sprintf('PSD of A: $\\sigma_A^2 = \\int S_{AA} df = %.3f$', var_A_freq), ...
    'Interpreter', 'latex')
xlim([0, 8])
grid on
box on
legend(band_labels, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 8)

% Panel 2: PSD of Signal B
subplot(2, 3, 2)
semilogy(f, Syy, 'Color', [0.8, 0.3, 0.2], 'LineWidth', 1.2)
hold on
for b = 1:num_bands
    idx = (f >= bands.(band_names{b})(1)) & (f < bands.(band_names{b})(2));
    area(f(idx), Syy(idx), 'FaceColor', band_colors{b}, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
end
hold off
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('$S_{BB}(f)$ [units$^2$/Hz]', 'Interpreter', 'latex')
title(sprintf('PSD of B: $\\sigma_B^2 = \\int S_{BB} df = %.3f$', var_B_freq), ...
    'Interpreter', 'latex')
xlim([0, 8])
grid on
box on

% Panel 3: Cross-spectral density (Real part)
subplot(2, 3, 3)
plot(f, real(Sxy), 'Color', [0.4, 0.2, 0.6], 'LineWidth', 1.2)
hold on
for b = 1:num_bands
    idx = (f >= bands.(band_names{b})(1)) & (f < bands.(band_names{b})(2));
    area(f(idx), real(Sxy(idx)), 'FaceColor', band_colors{b}, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
end
yline(0, '--', 'Color', [0.5, 0.5, 0.5])
hold off
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('Re$\{S_{AB}(f)\}$ [units$^2$/Hz]', 'Interpreter', 'latex')
title(sprintf('Cross-Spectral Density: Cov$(A,B) = \\int$ Re$\\{S_{AB}\\} df = %.3f$', cov_AB_freq), ...
    'Interpreter', 'latex')
xlim([0, 8])
grid on
box on

% Panel 4: Coherence
subplot(2, 3, 4)
plot(f, Cxy, 'k-', 'LineWidth', 1.2)
hold on
for b = 1:num_bands
    idx = (f >= bands.(band_names{b})(1)) & (f < bands.(band_names{b})(2));
    area(f(idx), Cxy(idx), 'FaceColor', band_colors{b}, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
end
hold off
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('$C_{AB}(f)$', 'Interpreter', 'latex')
title('Magnitude-Squared Coherence', 'Interpreter', 'latex')
xlim([0, 8])
% ylim([0, 1])
grid on
box on

% Panel 5: Phase spectrum
subplot(2, 3, 5)
plot(f, phi_xy * 180/pi, 'k-', 'LineWidth', 0.8)
hold on
for b = 1:num_bands
    idx = (f >= bands.(band_names{b})(1)) & (f < bands.(band_names{b})(2));
    % Highlight phase in bands where coherence is high
    high_coh = idx & (Cxy > 0.3);
    scatter(f(high_coh), phi_xy(high_coh)*180/pi, 20, 'black', 'filled')
end
yline(0, '--', 'Color', [0.5, 0.5, 0.5])
hold off
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('$\phi_{AB}(f)$ [deg]', 'Interpreter', 'latex')
title('Phase Spectrum (highlighted where $C_{AB} > 0.3$)', 'Interpreter', 'latex')
xlim([0, 8])
ylim([-180, 180])
yticks(-180:90:180)
grid on
box on

% Panel 6: Band contributions bar chart
subplot(2, 3, 6)
bar_data = zeros(num_bands, 3);
for b = 1:num_bands
    bar_data(b, :) = [band_results.(band_names{b}).frac_var_A, ...
                      band_results.(band_names{b}).frac_var_B, ...
                      band_results.(band_names{b}).frac_cov] * 100;
end
bh = bar(bar_data);
bh(1).FaceColor = [0.2, 0.4, 0.8];
bh(2).FaceColor = [0.8, 0.3, 0.2];
bh(3).FaceColor = [0.4, 0.2, 0.6];
set(gca, 'XTickLabel', {'LF', 'WF', 'HF'})
xlabel('Frequency Band', 'Interpreter', 'latex')
ylabel('Fraction of Total [\%]', 'Interpreter', 'latex')
title('Band Contributions to Variance \& Covariance', 'Interpreter', 'latex')
legend({'Var(A)', 'Var(B)', 'Cov(A,B)'}, 'Interpreter', 'latex', 'Location', 'northwest')
grid on
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3: PARSEVAL CONSERVATION DEMONSTRATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'white', 'Position', [150, 50, 1200, 700])
sgtitle('\textbf{Parseval''s Theorem: Band Decomposition Conserves Total Correlation}', ...
    'Interpreter', 'latex', 'FontSize', 16)

% Panel 1: Cumulative variance of A
subplot(2, 2, 1)
cum_Sxx = cumsum(Sxx) * df;
plot(f, cum_Sxx, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2)
hold on
yline(var_A, 'k--', 'LineWidth', 1.5, 'Label', sprintf('Var(A) = %.3f', var_A))
% Mark band boundaries
for b = 1:num_bands
    f_hi = bands.(band_names{b})(2);
    if f_hi < f(end)
        xline(f_hi, ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1)
    end
end
hold off
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('$\int_0^f S_{AA}(f'') df''$', 'Interpreter', 'latex')
title('Cumulative Integral of $S_{AA}$: Converges to Var(A)', 'Interpreter', 'latex')
xlim([0, fs/2])
grid on
box on

% Panel 2: Cumulative variance of B
subplot(2, 2, 2)
cum_Syy = cumsum(Syy) * df;
plot(f, cum_Syy, 'Color', [0.8, 0.3, 0.2], 'LineWidth', 2)
hold on
yline(var_B, 'k--', 'LineWidth', 1.5, 'Label', sprintf('Var(B) = %.3f', var_B))
for b = 1:num_bands
    f_hi = bands.(band_names{b})(2);
    if f_hi < f(end)
        xline(f_hi, ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1)
    end
end
hold off
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('$\int_0^f S_{BB}(f'') df''$', 'Interpreter', 'latex')
title('Cumulative Integral of $S_{BB}$: Converges to Var(B)', 'Interpreter', 'latex')
xlim([0, fs/2])
grid on
box on

% Panel 3: Cumulative covariance
subplot(2, 2, 3)
cum_Sxy_real = cumsum(real(Sxy)) * df;
plot(f, cum_Sxy_real, 'Color', [0.4, 0.2, 0.6], 'LineWidth', 2)
hold on
yline(cov_AB, 'k--', 'LineWidth', 1.5, 'Label', sprintf('Cov(A,B) = %.3f', cov_AB))
for b = 1:num_bands
    f_hi = bands.(band_names{b})(2);
    if f_hi < f(end)
        xline(f_hi, ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1)
    end
end
hold off
xlabel('Frequency [Hz]', 'Interpreter', 'latex')
ylabel('$\int_0^f \mathrm{Re}\{S_{AB}(f'')\} df''$', 'Interpreter', 'latex')
title('Cumulative Integral of Re$\{S_{AB}\}$: Converges to Cov(A,B)', 'Interpreter', 'latex')
xlim([0, fs/2])
grid on
box on

% Panel 4: Summary schematic
subplot(2, 2, 4)
axis off
hold on

% Draw boxes for each band
box_height = 0.15;
box_y = [0.7, 0.45, 0.2];
box_width = 0.25;

for b = 1:num_bands
    bname = band_names{b};
    
    % Band box
    rectangle('Position', [0.05, box_y(b), box_width, box_height], ...
        'FaceColor', band_colors{b}, 'EdgeColor', 'k', 'LineWidth', 1.5)
    
    % Band label
    text(0.05 + box_width/2, box_y(b) + box_height + 0.03, ...
        sprintf('%s [%.1f - %.1f Hz]', bname, bands.(bname)(1), bands.(bname)(2)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold')
    
    % Band statistics
    text(0.05 + box_width/2, box_y(b) + box_height/2, ...
        sprintf('$\\rho_{%s} = %.3f$', bname, band_results.(bname).rho), ...
        'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontSize', 12)
    
    % Arrow to sum
    annotation('arrow', [0.38, 0.48], [box_y(b) + box_height/2 + 0.12, 0.55], ...
        'LineWidth', 1.5, 'Color', band_colors{b})
end

% Sum box
rectangle('Position', [0.55, 0.35, 0.35, 0.3], 'FaceColor', [0.95, 0.95, 0.95], ...
    'EdgeColor', 'k', 'LineWidth', 2)
text(0.725, 0.72, '\textbf{Total Correlation}', 'HorizontalAlignment', 'center', ...
    'Interpreter', 'latex', 'FontSize', 12)
text(0.725, 0.55, sprintf('$\\rho = \\frac{\\sum_i \\mathrm{Cov}_i}{\\sqrt{\\sum_i \\mathrm{Var}_{A,i} \\cdot \\sum_i \\mathrm{Var}_{B,i}}}$'), ...
    'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontSize', 14)
text(0.725, 0.40, sprintf('$= %.4f$', rho_AB_freq), ...
    'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold')

% Title for panel
text(0.5, 0.95, '\textbf{Band Decomposition Framework}', ...
    'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontSize', 14)

hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4: BAND-SPECIFIC CORRELATION COEFFICIENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'white', 'Position', [200, 50, 1000, 400])
sgtitle('\textbf{Band-Partitioned Correlation Coefficients}', ...
    'Interpreter', 'latex', 'FontSize', 16)

% Bar chart of band correlations
subplot(1, 2, 1)
rho_bands = zeros(num_bands, 1);
for b = 1:num_bands
    rho_bands(b) = band_results.(band_names{b}).rho;
end

bh = bar(rho_bands, 'FaceColor', 'flat');
for b = 1:num_bands
    bh.CData(b, :) = sscanf(band_colors{b}(2:end), '%2x%2x%2x', [1 3]) / 255;
end
hold on
yline(rho_AB, 'k--', 'LineWidth', 2, 'Label', sprintf('Total ρ = %.3f', rho_AB))
yline(0, 'k-', 'LineWidth', 0.5)
hold off
set(gca, 'XTickLabel', band_labels)
ylabel('$\rho_{\mathrm{band}}$', 'Interpreter', 'latex')
title('Correlation Coefficient by Band', 'Interpreter', 'latex')
ylim([-0.5, 1])
grid on
box on

% Interpretation text
subplot(1, 2, 2)
axis off

text_str = {
    '\textbf{Interpretation:}', ...
    '', ...
    sprintf('$\\bullet$ \\textbf{LF band} ($\\rho = %.3f$): %s', band_results.LF.rho, ...
        interpret_rho(band_results.LF.rho)), ...
    '', ...
    sprintf('$\\bullet$ \\textbf{WF band} ($\\rho = %.3f$): %s', band_results.WF.rho, ...
        interpret_rho(band_results.WF.rho)), ...
    '', ...
    sprintf('$\\bullet$ \\textbf{HF band} ($\\rho = %.3f$): %s', band_results.HF.rho, ...
        interpret_rho(band_results.HF.rho)), ...
    '', ...
    '\textbf{Physical meaning:}', ...
    'High $\rho$ in WF band suggests turbines', ...
    'are synchronized by wave forcing.', ...
    'Low $\rho$ in LF band suggests independent', ...
    'slow drift motions.'
};

text(0.05, 0.9, text_str, 'Interpreter', 'latex', 'FontSize', 11, ...
    'VerticalAlignment', 'top', 'Units', 'normalized')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5: MATHEMATICAL FRAMEWORK SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'white', 'Position', [250, 50, 900, 700])

axis off
hold on

% Title
text(0.5, 0.95, '\textbf{Band-Partitioned Cross-Correlation: Mathematical Framework}', ...
    'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontSize', 16)

% Equations
eq_y = 0.85;
dy = 0.12;

text(0.05, eq_y, '\textbf{1. Parseval''s Theorem (Variance):}', ...
    'Interpreter', 'latex', 'FontSize', 12)
text(0.5, eq_y - 0.04, '$\displaystyle \sigma_x^2 = \mathrm{Var}(x) = \int_0^{f_s/2} S_{xx}(f) \, df = \sum_{\mathrm{bands}} \int_{\mathrm{band}} S_{xx}(f) \, df$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center')

text(0.05, eq_y - dy, '\textbf{2. Wiener-Khinchin Theorem (Covariance):}', ...
    'Interpreter', 'latex', 'FontSize', 12)
text(0.5, eq_y - dy - 0.04, '$\displaystyle \mathrm{Cov}(x,y) = \int_0^{f_s/2} \mathrm{Re}\{S_{xy}(f)\} \, df = \sum_{\mathrm{bands}} \int_{\mathrm{band}} \mathrm{Re}\{S_{xy}(f)\} \, df$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center')

text(0.05, eq_y - 2*dy, '\textbf{3. Total Correlation Coefficient:}', ...
    'Interpreter', 'latex', 'FontSize', 12)
text(0.5, eq_y - 2*dy - 0.04, '$\displaystyle \rho_{xy} = \frac{\mathrm{Cov}(x,y)}{\sqrt{\mathrm{Var}(x) \cdot \mathrm{Var}(y)}} = \frac{\int \mathrm{Re}\{S_{xy}\} df}{\sqrt{\int S_{xx} df \cdot \int S_{yy} df}}$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center')

text(0.05, eq_y - 3*dy, '\textbf{4. Band-Partitioned Correlation:}', ...
    'Interpreter', 'latex', 'FontSize', 12)
text(0.5, eq_y - 3*dy - 0.04, '$\displaystyle \rho_{\mathrm{band}} = \frac{\int_{\mathrm{band}} \mathrm{Re}\{S_{xy}(f)\} \, df}{\sqrt{\int_{\mathrm{band}} S_{xx}(f) \, df \cdot \int_{\mathrm{band}} S_{yy}(f) \, df}}$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center')

text(0.05, eq_y - 4*dy, '\textbf{5. Key Property (Conservation):}', ...
    'Interpreter', 'latex', 'FontSize', 12)
text(0.5, eq_y - 4*dy - 0.04, '$\displaystyle \sum_{\mathrm{bands}} \mathrm{Var}_{\mathrm{band}} = \mathrm{Var}_{\mathrm{total}}, \quad \sum_{\mathrm{bands}} \mathrm{Cov}_{\mathrm{band}} = \mathrm{Cov}_{\mathrm{total}}$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center')

% Note box
rectangle('Position', [0.05, 0.02, 0.9, 0.15], 'FaceColor', [0.95, 0.95, 1], ...
    'EdgeColor', [0.2, 0.2, 0.6], 'LineWidth', 1.5, 'Curvature', 0.1)
text(0.5, 0.12, '\textbf{Note:} $\rho_{\mathrm{band}}$ values do NOT sum to $\rho_{\mathrm{total}}$ because correlation is a normalized quantity.', ...
    'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontSize', 11)
text(0.5, 0.06, 'However, the \textit{unnormalized} quantities (variance, covariance) DO sum, satisfying Parseval''s theorem.', ...
    'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontSize', 11)

hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINT FINAL SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
fprintf('=========================================================\n')
fprintf('   FINAL SUMMARY: PARSEVAL CONSERVATION CHECK\n')
fprintf('=========================================================\n\n')

sum_var_A = 0; sum_var_B = 0; sum_cov = 0;
for b = 1:num_bands
    sum_var_A = sum_var_A + band_results.(band_names{b}).var_A;
    sum_var_B = sum_var_B + band_results.(band_names{b}).var_B;
    sum_cov = sum_cov + band_results.(band_names{b}).cov_AB;
end

fprintf('Sum of band Var(A):   %.6f  (Total: %.6f, Error: %.4f%%)\n', ...
    sum_var_A, var_A_freq, 100*abs(sum_var_A - var_A_freq)/var_A_freq)
fprintf('Sum of band Var(B):   %.6f  (Total: %.6f, Error: %.4f%%)\n', ...
    sum_var_B, var_B_freq, 100*abs(sum_var_B - var_B_freq)/var_B_freq)
fprintf('Sum of band Cov(A,B): %.6f  (Total: %.6f, Error: %.4f%%)\n', ...
    sum_cov, cov_AB_freq, 100*abs(sum_cov - cov_AB_freq)/abs(cov_AB_freq))
fprintf('\n')
fprintf('Reconstructed ρ from band sums: %.6f\n', sum_cov / sqrt(sum_var_A * sum_var_B))
fprintf('Direct frequency-domain ρ:      %.6f\n', rho_AB_freq)
fprintf('Time-domain ρ:                  %.6f\n', rho_AB)
fprintf('\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = interpret_rho(rho)
    if rho > 0.7
        desc = 'Strong positive coupling';
    elseif rho > 0.3
        desc = 'Moderate positive coupling';
    elseif rho > -0.3
        desc = 'Weak/negligible coupling';
    elseif rho > -0.7
        desc = 'Moderate negative coupling';
    else
        desc = 'Strong negative coupling';
    end
end