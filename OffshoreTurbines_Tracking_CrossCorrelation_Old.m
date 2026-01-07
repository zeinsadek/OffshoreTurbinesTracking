%%% Offshore Tracking Signals

clear; close all; clc;
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT TRACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

offshore_path = "/Users/zeinsadek/Desktop/Experiments/Offshore";
power_path = fullfile(offshore_path, "Power/Data/Matfiles");
tracking_path = fullfile(offshore_path, "Tracking/Data/Matfiles");

farm_arrangement = "Inline";
farm_spacing = "SX50";
caze = strcat("WT60_", farm_spacing, "_AG0");

% Import Tracking
tracking_file = fullfile(tracking_path, farm_arrangement, strcat(caze, ".mat"));
tmp_tracking = load(tracking_file);
tmp_tracking = tmp_tracking.(caze);

clear tracking_path tracking_file

% Convert fieldnames into arry for turbines in tracking structure to match
% power
waves = fieldnames(tmp_tracking);
turbines = fieldnames(tmp_tracking.(waves{1}));

for w = 1:length(waves)
    wave = waves{w};
    for t = 1:length(turbines)
        turbine = turbines{t};
        tracking.(wave)(t) = tmp_tracking.(wave).(turbine);
    end
end

clear wave turbine w t tmp_tracking


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVE FORCING FREQUENCIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

for w = 1:length(waves)
    wave = waves{w};
    for t = 1:length(turbines)

        % Time
        corrected_tracking.(wave)(t).time = tracking.(wave)(t).time;

        % Raw
        corrected_tracking.(wave)(t).x = fillmissing(tracking.(wave)(t).x, fillmethod);
        corrected_tracking.(wave)(t).y = fillmissing(tracking.(wave)(t).z, fillmethod);
        corrected_tracking.(wave)(t).z = fillmissing(-1 * tracking.(wave)(t).y, fillmethod);

        corrected_tracking.(wave)(t).roll = fillmissing(tracking.(wave)(t).roll, fillmethod);
        corrected_tracking.(wave)(t).pitch = fillmissing(-1 * tracking.(wave)(t).pitch, fillmethod);
        corrected_tracking.(wave)(t).yaw = fillmissing(tracking.(wave)(t).yaw, fillmethod);
        
        % Kalman
        corrected_tracking.(wave)(t).x_kal = fillmissing(tracking.(wave)(t).x_kal, fillmethod);
        corrected_tracking.(wave)(t).y_kal = fillmissing(tracking.(wave)(t).z_kal, fillmethod);
        corrected_tracking.(wave)(t).z_kal = fillmissing(-1 * tracking.(wave)(t).y_kal, fillmethod);

        corrected_tracking.(wave)(t).roll_kal = fillmissing(tracking.(wave)(t).roll_kal, fillmethod);
        corrected_tracking.(wave)(t).pitch_kal = fillmissing(-1 * tracking.(wave)(t).pitch_kal, fillmethod);
        corrected_tracking.(wave)(t).yaw_kal = fillmissing(tracking.(wave)(t).yaw_kal, fillmethod);

    end
end


% Replace with correct signs
clear tracking w t wave
tracking = corrected_tracking;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sliding mean window size
window = 6;

% Colors
colors = {'#710627', '#EF6461', '#9CE37D', '#8ACDEA'};

% Plot details
linewidth = 1.5;
marker_size = 15;

% Which wave
wave = 'LM4_AK12';
wave_parts = split(wave, '_');
wavelength = wave_parts{1};

% Degree of freedom
center_turbines = [2,5,8,10];
DOF = 'pitch_kal';

% Plot
clc; close all;
ax = figure('color', 'white');
tiledlayout(4,1)
sgtitle(wave, 'interpreter', 'none')
for i = 1:length(center_turbines)
    nexttile
    turbine = center_turbines(i);

    hold on
    plot(tracking.(wave)(turbine).time, movmean(tracking.(wave)(turbine).(DOF), window), ...
         'Color', colors{i}, 'LineWidth', linewidth)
    
    scatter(tracking.(wave)(turbine).time, movmean(tracking.(wave)(turbine).(DOF), window), ...
            marker_size, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors{i})
    hold off

    title(strcat("Row", " ", num2str(i)), 'Interpreter', 'None')
    xlim([10, 20])
    ylim([-30, 10])

    if i == 2
        ylabel('Pitch [Degrees]')
    end

    if i == 4
        xlabel('Time [Seconds]')
    end
    clear i
end

% Save fig
% exportgraphics(ax, fullfile(out_path, strcat(wave, '_Pitch_Sample.png')), 'Resolution', 200)


%%% Plot w signals overlayed
ax = figure('color', 'white');
sgtitle(wave, 'interpreter', 'none')
hold on
for i = 1:length(center_turbines)
    turbine = center_turbines(i);
    plot(tracking.(wave)(turbine).time, movmean(tracking.(wave)(turbine).(DOF), window), ...
         'Color', colors{i}, 'LineWidth', linewidth, 'DisplayName', strcat('Row', " ", num2str(i)))
    
    scatter(tracking.(wave)(turbine).time, movmean(tracking.(wave)(turbine).(DOF), window), ...
            marker_size, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors{i}, ...
            'HandleVisibility', 'off') 
end
hold off
xlim([16, 17])  
ylim([-30, 10])
legend()
ylabel('Pitch [Degrees]')
xlabel('Time [Seconds]')

% Save fig
% exportgraphics(ax, fullfile(out_path, strcat(wave, '_Pitch_Overlayed.png')), 'Resolution', 200)


%% % Cross Correlation

% Time steps
dt = mean(diff(tracking.(wave)(1).time), 'omitnan');    

% Wave period
wave_period = 1 / forcing_frequencies.(wavelength);

% Crop signal
end_buffer = 300;

% Harmonic
if strcmp(wavelength, 'LM5')
    H = 1;
elseif strcmp(wavelength, 'LM4')
    H = 1.25;
elseif strcmp(wavelength, 'LM33')
    H = 1.5;
elseif strcmp(wavelength, 'LM25')
    H = 2;
end

% Which degree of freedom
DOF = 'pitch_kal';

% Plot
ax = figure('Position', [200,200,800,400], 'color', 'white');

hold on
for i = 1:length(center_turbines) - 1
    % Front turbine
    signal1 = movmean(tracking.(wave)(center_turbines(1)).(DOF)(1:end - end_buffer), window);
    signal1 = signal1 - mean(signal1);

    % Waked turbines
    signal2 = movmean(tracking.(wave)(center_turbines(i + 1)).(DOF)(1:end - end_buffer), window);
    signal2 = signal2 - mean(signal2);

    % Cross-correlation
    [cor, lags] = xcorr(signal1, signal2, 'normalized');

    % Convert steps into time
    lags = (lags * dt) / wave_period;

    % Plot
    plot(lags, cor, 'DisplayName', strcat('Row 1 X Row', " " , num2str(i + 1)), ...
         'linewidth', 3, 'Color', colors{i + 1})

    clear i signal1 signal2 corr lags
end
hold on
ylim([-1, 1])
xlim([0, 10])
title(wave, 'Interpreter','none')
xlabel('$t / T_{wave}$', 'interpreter', 'latex')
ylabel('Normalized Correlation')
legend()
% exportgraphics(ax, fullfile(out_path, strcat(wave, '_Pitch_Cross_Correlation.png')), 'Resolution', 200)


%% Loop through different waves and generate a tiledlayout



waves = {'LM5_AK12', 'LM4_AK12', 'LM33_AK12', 'LM25_AK12'};
wave_labels = {"H = 1", "H = 1.25", "H = 1.5", "H = 2"};



% Time steps
dt = mean(diff(tracking.(wave)(1).time), 'omitnan');    

% Wave period
wave_period = 1 / forcing_frequencies.(wavelength);

% Crop signal
end_buffer = 300;

% Harmonic
if strcmp(wavelength, 'LM5')
    H = 1;
elseif strcmp(wavelength, 'LM4')
    H = 1.25;
elseif strcmp(wavelength, 'LM33')
    H = 1.5;
elseif strcmp(wavelength, 'LM25')
    H = 2;
end

% Which degree of freedom
DOF = 'pitch_kal';

% Plot
ax = figure('Position', [200,200,800,400], 'color', 'white');
tiledlayout(length(waves), 1)

for w = 1:length(waves)
    h(w) = nexttile;
    wave = waves{w};

    hold on
    for i = 1:length(center_turbines) - 1
        % Front turbine
        signal1 = movmean(tracking.(wave)(center_turbines(1)).(DOF)(1:end - end_buffer), window);
        signal1 = signal1 - mean(signal1);
    
        % Waked turbines
        signal2 = movmean(tracking.(wave)(center_turbines(i + 1)).(DOF)(1:end - end_buffer), window);
        signal2 = signal2 - mean(signal2);
    
        % Cross-correlation
        [cor, lags] = xcorr(signal1, signal2, 'normalized');
    
        % Convert steps into time
        lags = (lags * dt) / wave_period;
    
        % Plot
        plot(lags, cor, 'DisplayName', strcat('Row 1 X Row', " " , num2str(i + 1)), ...
             'linewidth', 3, 'Color', colors{i + 1})
    
        clear i signal1 signal2 corr lags
    end
    hold off
    ylim([-1, 1])
    xlim([0, 1])
    title(wave_labels{w}, 'Interpreter','latex')
    ylabel('Normalized Correlation')
    if w == 4
        xlabel('$t / T_{wave}$', 'interpreter', 'latex')
    end
    if w == 1
        legend('location', 'northeastoutside')
    end
end

linkaxes(h, 'xy')