%% Visualization of Band-Partitioned Variance Decomposition
% Zein Sadek
%
% This script creates a pedagogical figure showing:
% 1. The original time-domain signal
% 2. The power spectral density with shaded frequency bands
% 3. A bar chart showing variance partitioning and Parseval's theorem
%
% Intended for inclusion in papers/presentations to explain the methodology

clear; close all; clc;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT TRACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

offshore_path = "/Users/zeinsadek/Desktop/Experiments/Offshore";
power_path = fullfile(offshore_path, "Power/Data/Matfiles");
tracking_path = fullfile(offshore_path, "Tracking/Data/Matfiles");

farm_arrangement = "Inline";

% Load all spacings
farm_spacings = [5, 4.5, 4, 3.5, 3];

for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
    caze = strcat("WT60_", farm_spacing, "_AG0");
    fprintf('%s\n', caze)
    
    % Import Tracking
    tracking_file = fullfile(tracking_path, farm_arrangement, strcat(caze, ".mat"));
    tmp = load(tracking_file);
    tmp = tmp.(caze);

    tmp_tracking.(farm_spacing) = tmp;

    clear s caze tracking_file tmp
end
clear tracking_path 

% Convert fieldnames into arry for turbines in tracking structure to match
% power
for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
    caze = strcat("WT60_", farm_spacing, "_AG0");
    fprintf('%s\n', caze)
    tmp = tmp_tracking.(farm_spacing);

    waves = fieldnames(tmp);
    turbines = fieldnames(tmp.(waves{1}));

    for w = 1:length(waves)
        wave = waves{w};
        disp(wave)
        for t = 1:length(turbines)
            turbine = turbines{t};
            tracking.(farm_spacing).(wave)(t) = tmp_tracking.(farm_spacing).(wave).(turbine);
        end

    end
end

clear wave turbine w t tmp_tracking

forcing_frequencies.("LM5")  = 1.4273;
forcing_frequencies.("LM4")  = 1.6075;
forcing_frequencies.("LM33") = 1.7651;
forcing_frequencies.("LM3")  = 1.8617;
forcing_frequencies.("LM25") = 2.0402;
forcing_frequencies.("LM2")  = 2.2813;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECT COORDINATE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fillmethod = 'spline';

clc;
for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
    caze = strcat("WT60_", farm_spacing, "_AG0");
    fprintf('%s\n', caze)
    
    tmp = tracking.(farm_spacing);
    waves = fieldnames(tmp);

    for w = 1:length(waves)
        wave = waves{w};
        for t = 1:length(turbines)
    
            % Time
            corrected_tracking.(farm_spacing).(wave)(t).time = tracking.(farm_spacing).(wave)(t).time;
    
            % Raw
            % Convert cm to m
            corrected_tracking.(farm_spacing).(wave)(t).x = fillmissing(tracking.(farm_spacing).(wave)(t).x, fillmethod) .* 1E-2;
            corrected_tracking.(farm_spacing).(wave)(t).y = fillmissing(tracking.(farm_spacing).(wave)(t).z, fillmethod) .* 1E-2;
            corrected_tracking.(farm_spacing).(wave)(t).z = fillmissing(-1 * tracking.(farm_spacing).(wave)(t).y, fillmethod) .* 1E-2;
    
            % Rotation in degrees
            corrected_tracking.(farm_spacing).(wave)(t).roll = fillmissing(tracking.(farm_spacing).(wave)(t).roll, fillmethod);
            corrected_tracking.(farm_spacing).(wave)(t).pitch = fillmissing(-1 * tracking.(farm_spacing).(wave)(t).pitch, fillmethod);
            corrected_tracking.(farm_spacing).(wave)(t).yaw = fillmissing(tracking.(farm_spacing).(wave)(t).yaw, fillmethod);
            
            % Kalman
            % Convert cm to m
            corrected_tracking.(farm_spacing).(wave)(t).x_kal = fillmissing(tracking.(farm_spacing).(wave)(t).x_kal, fillmethod) .* 1E-2;
            corrected_tracking.(farm_spacing).(wave)(t).y_kal = fillmissing(tracking.(farm_spacing).(wave)(t).z_kal, fillmethod) .* 1E-2;
            corrected_tracking.(farm_spacing).(wave)(t).z_kal = fillmissing(-1 * tracking.(farm_spacing).(wave)(t).y_kal, fillmethod) .* 1E-2;
    
            % Rotation in degrees
            corrected_tracking.(farm_spacing).(wave)(t).roll_kal = fillmissing(tracking.(farm_spacing).(wave)(t).roll_kal, fillmethod);
            corrected_tracking.(farm_spacing).(wave)(t).pitch_kal = fillmissing(-1 * tracking.(farm_spacing).(wave)(t).pitch_kal, fillmethod);
            corrected_tracking.(farm_spacing).(wave)(t).yaw_kal = fillmissing(tracking.(farm_spacing).(wave)(t).yaw_kal, fillmethod);
    
        end
    end
end


% Replace with correct signs
clear tracking w t
tracking = corrected_tracking;



%% ==========================================================
% USER PARAMETERS
% ===========================================================

% Figure sizing
figWidth  = 7.5;   % inches (fits journal double-column)
figHeight = 7;     % inches

% Font sizes
titleFontSize  = 11;
labelFontSize  = 10;
tickFontSize   = 9;
annotFontSize  = 9;

% Export settings
exportFig  = false;
exportPath = 'BandPartitionedVariance_Schematic.pdf';

% Signal parameters (adjust to match your data)
Fs = 30;           % Sampling frequency [Hz]
fw = 1.4273;       % Wave forcing frequency [Hz] (e.g., LM5)

% Band definitions
bands.LF   = [0,       0.5*fw];
bands.Wave = [0.5*fw,  1.5*fw];
bands.HF   = [1.5*fw,  Fs/2];

% Colors for frequency bands
colors.LF   = [0.400, 0.761, 0.647];  % Teal/green - low frequency
colors.Wave = [0.988, 0.553, 0.384];  % Coral/orange - wave frequency
colors.HF   = [0.553, 0.627, 0.796];  % Steel blue - high frequency
colors.total = [0.3, 0.3, 0.3];       % Dark gray for total



%% ==========================================================
% LOAD OR GENERATE EXAMPLE SIGNAL
% ===========================================================

% Option 1: Load your actual data
% Uncomment and modify the path to use real data:
% -------------------------
% load('your_tracking_data.mat');
% t = tracking.SX50.LM5_AK12(1).time;
% q = tracking.SX50.LM5_AK12(1).pitch_kal;
% -------------------------

% Option 2: Generate synthetic signal with known frequency content
% This creates a signal with clear LF, Wave, and HF components
% T = 60;                         % Duration [s]
% t = (0:1/Fs:T-1/Fs)';           % Time vector
% N = length(t);

% Test cross-correlation
Fs = 30;
farm_spacing = 'SX50';
wave = 'LM5_AK12';
fw = forcing_frequencies.LM5;
DOF = 'pitch_kal';

% Signal
t = tracking.(farm_spacing).(wave)(2).time;
N = length(t);
q = tracking.(farm_spacing).(wave)(2).(DOF);


% Remove mean (as we do in analysis)
q = q - mean(q);

% fprintf('Signal: N = %d samples, T = %.1f s, Fs = %.1f Hz\n', N, T, Fs);
% fprintf('Wave frequency: fw = %.4f Hz\n', fw);


%% ==========================================================
% COMPUTE PSD AND BAND-LIMITED VARIANCES
% ===========================================================

% Compute FFT
Q = fft(q);
df = Fs / N;
f = (0:N-1)' * df;

% One-sided frequency vector
if mod(N, 2) == 0
    nUnique = N/2 + 1;
else
    nUnique = (N+1)/2;
end
f_onesided = f(1:nUnique);
Q_onesided = Q(1:nUnique);

% Compute one-sided PSD (power per Hz)
Pxx = abs(Q_onesided).^2 / (N^2 * df);
% Scale for one-sided (×2 except DC and Nyquist)
Pxx(2:end-1) = 2 * Pxx(2:end-1);
if mod(N, 2) ~= 0
    Pxx(end) = 2 * Pxx(end);
end

% Find indices for each band
idx_LF   = (f_onesided >= bands.LF(1))   & (f_onesided < bands.LF(2));
idx_Wave = (f_onesided >= bands.Wave(1)) & (f_onesided < bands.Wave(2));
idx_HF   = (f_onesided >= bands.HF(1))   & (f_onesided <= bands.HF(2));

% Compute variance in each band (integral of PSD)
var_LF   = trapz(f_onesided(idx_LF), Pxx(idx_LF));
var_Wave = trapz(f_onesided(idx_Wave), Pxx(idx_Wave));
var_HF   = trapz(f_onesided(idx_HF), Pxx(idx_HF));

% Total variance from signal
var_total = var(q);

% Sum of band variances
var_sum = var_LF + var_Wave + var_HF;

% RMS values
rms_LF   = sqrt(var_LF);
rms_Wave = sqrt(var_Wave);
rms_HF   = sqrt(var_HF);
rms_total = std(q);

fprintf('\n--- Variance Decomposition ---\n');
fprintf('LF variance:   %.4f  (RMS = %.4f)\n', var_LF, rms_LF);
fprintf('Wave variance: %.4f  (RMS = %.4f)\n', var_Wave, rms_Wave);
fprintf('HF variance:   %.4f  (RMS = %.4f)\n', var_HF, rms_HF);
fprintf('Sum of bands:  %.4f\n', var_sum);
fprintf('Total variance: %.4f  (STD = %.4f)\n', var_total, rms_total);
fprintf('Closure ratio: %.6f\n', var_sum / var_total);


%% ==========================================================
% CREATE FIGURE
% ===========================================================

% Layout: 3 rows
% Row 1: Time series
% Row 2: PSD with shaded bands
% Row 3: Variance bar chart

%% ----- Panel (a): Time Series -----


figure('color', 'white')
hold on;
plot(t, q, 'k-', 'LineWidth', 2);

xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('$q(t)$', 'Interpreter', 'latex', 'FontSize', labelFontSize);
title('(a) Example Signal: Pitch Motion', 'Interpreter', 'latex', 'FontSize', titleFontSize);

% Add variance annotation
text(0.98, 0.92, sprintf('$\\sigma_q = %.2f$', rms_total), ...
    'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', annotFontSize, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);

xlim([0, min(20, 30)]);  % Show first 20 seconds for clarity
set(gca, 'FontSize', tickFontSize, 'TickLabelInterpreter', 'latex');
box on; grid on;
hold off;


%% ----- Panel (b): Power Spectral Density with Shaded Bands -----

figure('color', 'white')
hold on;
% Shade frequency bands (must be done before plotting PSD for proper layering)
ymax_shade = max(Pxx) * 1.1;

% LF band
fill([bands.LF(1)/fw, bands.LF(2)/fw, bands.LF(2)/fw, bands.LF(1)/fw], ...
     [0, 0, ymax_shade, ymax_shade], ...
     colors.LF, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Wave band
fill([bands.Wave(1)/fw, bands.Wave(2)/fw, bands.Wave(2)/fw, bands.Wave(1)/fw], ...
     [0, 0, ymax_shade, ymax_shade], ...
     colors.Wave, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% HF band
fill([bands.HF(1)/fw, bands.HF(2)/fw, bands.HF(2)/fw, bands.HF(1)/fw], ...
     [0, 0, ymax_shade, ymax_shade], ...
     colors.HF, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot PSD
plot(f_onesided/fw, Pxx, 'k-', 'LineWidth', 1);

% Mark band boundaries with vertical lines
xline(bands.LF(2)/fw, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
xline(bands.Wave(2)/fw, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);

% Mark wave frequency
% xline(fw, ':', 'Color', colors.Wave*0.7, 'LineWidth', 1.5);
text(fw, ymax_shade*0.95, '$f_w$', 'Interpreter', 'latex', ...
    'FontSize', annotFontSize, 'HorizontalAlignment', 'center', ...
    'Color', colors.Wave*0.6);

% Add band labels
text(mean(bands.LF), ymax_shade*0.85, 'LF', 'FontSize', annotFontSize, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', colors.LF*0.7);
text(mean(bands.Wave), ymax_shade*0.85, 'Wave', 'FontSize', annotFontSize, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', colors.Wave*0.7);
text(mean([bands.HF(1), min(bands.HF(2), 8)]), ymax_shade*0.85, 'HF', 'FontSize', annotFontSize, ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', colors.HF*0.7);

% Add integral annotations
text(mean(bands.LF), ymax_shade*0.72, sprintf('$\\sigma_{LF}^2 = %.2f$', var_LF), ...
    'Interpreter', 'latex', 'FontSize', annotFontSize-1, ...
    'HorizontalAlignment', 'center', 'Color', colors.LF*0.6);
text(mean(bands.Wave), ymax_shade*0.72, sprintf('$\\sigma_{W}^2 = %.2f$', var_Wave), ...
    'Interpreter', 'latex', 'FontSize', annotFontSize-1, ...
    'HorizontalAlignment', 'center', 'Color', colors.Wave*0.6);

xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel('PSD [$\cdot^2$/Hz]', 'Interpreter', 'latex', 'FontSize', labelFontSize);
title('(b) Power Spectral Density with Frequency Band Partitioning', ...
    'Interpreter', 'latex', 'FontSize', titleFontSize);

xlim([0, min(Fs/2, 5)]);  % Limit x-axis for clarity
ylim([0, ymax_shade]);
set(gca, 'FontSize', tickFontSize, 'TickLabelInterpreter', 'latex');
box on;
hold off;
% yscale('log')
% ylim([0, 1E-3])

%% ----- Panel (c): Variance Bar Chart -----


figure('color', 'white')
hold on;

% Bar positions
barX = [1, 2, 3, 4.5, 6];
barLabels = {'LF', 'Wave', 'HF', '$\Sigma$', 'Total'};

% Bar values (variances)
barY = [var_LF, var_Wave, var_HF, var_sum, var_total];

% Bar colors
barColors = [colors.LF; colors.Wave; colors.HF; [0.5 0.5 0.5]; colors.total];

% Draw bars
for i = 1:length(barX)
    bar(barX(i), barY(i), 0.6, 'FaceColor', barColors(i,:), 'EdgeColor', 'k', 'LineWidth', 1);
end

% Add "+" signs between first three bars
text(1.5, max(barY)*0.15, '+', 'FontSize', 14, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(2.5, max(barY)*0.15, '+', 'FontSize', 14, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Add "=" sign before sum
text(3.75, max(barY)*0.15, '=', 'FontSize', 14, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Add "≈" between sum and total
text(5.25, max(barY)*0.15, '$\approx$', 'FontSize', 14, 'HorizontalAlignment', 'center', ...
    'Interpreter', 'latex', 'FontWeight', 'bold');

% Add value labels on bars
for i = 1:length(barX)
    text(barX(i), barY(i) + max(barY)*0.03, sprintf('%.2f', barY(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', annotFontSize-1);
end

% Add closure ratio annotation
text(5.5, max(barY)*0.85, sprintf('Closure: %.4f', var_sum/var_total), ...
    'FontSize', annotFontSize, 'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);

% Add equation
text(0.5, 0.95, '$\sigma_{LF}^2 + \sigma_{W}^2 + \sigma_{HF}^2 = \sigma_{total}^2$ \quad (Parseval''s Theorem)', ...
    'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', annotFontSize, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

xlabel('', 'FontSize', labelFontSize);
ylabel('Variance [$\cdot^2$]', 'Interpreter', 'latex', 'FontSize', labelFontSize);
title('(c) Variance Partitioning Verification', 'Interpreter', 'latex', 'FontSize', titleFontSize);

set(gca, 'XTick', barX, 'XTickLabel', barLabels, 'TickLabelInterpreter', 'latex');
xlim([0.3, 6.7]);
ylim([0, max(barY)*1.15]);
set(gca, 'FontSize', tickFontSize);
box on;
hold off;


% %% ==========================================================
% % EXPORT FIGURE
% % ===========================================================
% 
% if exportFig
%     set(fig, 'PaperUnits', 'inches');
%     set(fig, 'PaperSize', [figWidth, figHeight]);
%     set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
% 
%     if endsWith(exportPath, '.pdf')
%         exportgraphics(fig, exportPath, 'ContentType', 'vector');
%     elseif endsWith(exportPath, '.png')
%         exportgraphics(fig, exportPath, 'Resolution', 300);
%     elseif endsWith(exportPath, '.eps')
%         exportgraphics(fig, exportPath, 'ContentType', 'vector');
%     end
% 
%     fprintf('\nFigure exported to: %s\n', exportPath);
% end
% 
% fprintf('\n=== Figure Generation Complete ===\n');