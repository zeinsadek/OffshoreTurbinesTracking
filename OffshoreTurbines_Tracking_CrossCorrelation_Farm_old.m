%% Looking at cross-correlation of motion between turbines/rows

% Cross-correlating different rows of turbines, forcing positive-definite
% time lags since waked rows have to respond to first row, they cannot
% anticipate its motion

% Zein Sadek

clear; close all; clc;
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM")
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT TRACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Paths
offshore_path = "/Users/zeinsadek/Desktop/Experiments/Offshore";
tracking_path = fullfile(offshore_path, "Tracking/Data/Matfiles");

% Farm layout + Data filter version
farm_arrangement = "Staggered";
harmonic_cutoff = 2;


% Number of turbines based on arrangement
turbine_catalog.Inline.turbines = 1:12;
turbine_catalog.Staggered.turbines = 1:10;

% Center turbines based on arrangement
turbine_catalog.Inline.centers = [2,5,8,10];
turbine_catalog.Staggered.centers = [2,4,7,9];

% Get farm layout info
turbines = turbine_catalog.(farm_arrangement).turbines;
centers = turbine_catalog.(farm_arrangement).centers;

% Farm spacings
farm_spacings = [5, 4.5, 4, 3.5, 3];

% Load tracking
tracking_file = sprintf('OffshoreTracking_AllDataCombined_SavitskyGolay_Cutoff_%1.0f.mat', harmonic_cutoff);
full_tracking = load(fullfile(tracking_path, "AllData", tracking_file));
full_tracking = full_tracking.tracking;

% Load Savitsky-Golay filter design
filterspecs = full_tracking.FilterDesign;

% Load data for specific layout
tracking = full_tracking.(farm_arrangement);

clear offshore_path tracking_file


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVE INFORMATION AND NAMING CONVENTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Farm arrangements
farm_arrangements = {"Inline", "Staggered"};

% Wave forcing frequencies
forcing_frequencies.("LM5")  = 1.4273;
forcing_frequencies.("LM4")  = 1.6075;
forcing_frequencies.("LM33") = 1.7651;
forcing_frequencies.("LM3")  = 1.8617;
forcing_frequencies.("LM25") = 2.0402;
forcing_frequencies.("LM2")  = 2.2813;
% Wavelengths
wavelengths_from_string.("LM5")  = 5;
wavelengths_from_string.("LM4")  = 4;
wavelengths_from_string.("LM33") = (5 / 1.5);
wavelengths_from_string.("LM3")  = 3;
wavelengths_from_string.("LM25") = 2.5;
wavelengths_from_string.("LM2")  = 2;
% Steepnesses
steepesses_from_string.("AK12") = 0.12;
steepesses_from_string.("AK09") = 0.09;
steepesses_from_string.("AK06") = 0.06;


% Degrees of freedom to loop over
DOFs = {'x_kal', 'y_kal', 'z_kal', 'roll_kal', 'pitch_kal', 'yaw_kal'};
% Plotting nomenclature: DOF names
DOF_names.x_kal = 'Surge';
DOF_names.y_kal = 'Heave';
DOF_names.z_kal = 'Sway';
DOF_names.roll_kal = 'Roll';
DOF_names.pitch_kal = 'Pitch';
DOF_names.yaw_kal = 'Yaw';
% Plotting nomenclature: DOF symbols
DOF_symbs.x_kal = 'x';
DOF_symbs.y_kal = 'y';
DOF_symbs.z_kal = 'z';
DOF_symbs.roll_kal = '\phi';
DOF_symbs.pitch_kal = '\theta';
DOF_symbs.yaw_kal = '\psi';
% Types of motion
translations = {'x_kal', 'y_kal', 'z_kal'};
rotations = {'roll_kal', 'pitch_kal', 'yaw_kal'};

% Experiment details
sampling_frequency = 30;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING PROPERTIES AND CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wavelengths and steepnesses to loop over
wavelengths = [5,4,3,2];
wave_steepnesses = [0.06, 0.09, 0.12];

% Colors per row
row_colors.Row1 = flipud(slanCM(45, 2 * length(wavelengths)));
row_colors.Row2 = flipud(slanCM(47, 2 * length(wavelengths)));
row_colors.Row3 = flipud(slanCM(34, 2 * length(wavelengths)));
row_colors.Row4 = flipud(slanCM(35, 2 * length(wavelengths)));

% Scatter plot nomenclature
marker_size = 100;
steepness_alpha = [0.3, 0.6, 1];
spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ALL CROSS-CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG: MAX CORRELATION VALUE (WITH SIGN) AGAINST
% HARMOIC RATIO
% LOOPED OVER WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correlating a specified DOF from a fixed row, to all waked rows
% Plotting the absolute value of maximum correlation coefficient

% Which row to correlate against
fixed_row = 1;
DOF = 'yaw_kal';
colors = row_colors.(sprintf('Row%1.0f', fixed_row));


% Cross-correlation 
wave_advection_multiplier = 1.5;

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);    

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plotting
clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('Maximum Correlation Coefficient of %s ($%s$): $|\\rho_{%s}|$', name, symb, symb), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    % Get turbines
    reference_turbine = centers(fixed_row);
    correlating_turbine = waked_turbines(c);

    h(c) = nexttile;
    title(sprintf('Row %1.0f to Row %1.0f', fixed_row, ceil(correlating_turbine / 3)))
    hold on
    % Loop through spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through waves
        for st = 1:length(wave_steepnesses)
            steep = compose('%02d', round(100 * wave_steepnesses(st)));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)

                    % How many waves fit between rows
                    harmonic_ratio = farm_spacings(s) / wavelengths(w);

                    % How long it takes for wave to travel between rows
                    phase_velocity = wavelengths(w) * forcing_frequencies.(['LM', num2str(wavelengths(w))]);
                    wave_advection_time = farm_spacings(s) / phase_velocity;
                    maxLagSeconds = wave_advection_multiplier * wave_advection_time;
    
                    % Get signals
                    turbine_A_signal = tracking.(farm_spacing).(wave)(reference_turbine).(DOF);
                    turbine_B_signal = tracking.(farm_spacing).(wave)(correlating_turbine).(DOF);
                    
                    % Make the same length since we use 'coeff'
                    n = min(numel(turbine_A_signal), numel(turbine_B_signal));
                    A = turbine_A_signal(1:n);
                    B = turbine_B_signal(1:n);
        
                    % Cross-correlate
                    out = xcorr_metrics(A, B, sampling_frequency, maxLagSeconds);
        
                    % out.rho_max ~ largest magnitude XC coefficient
                    % out.rho_abs_max ~ abs of largest magnitude XC coefficient
                    % out.tau_max ~ time lag at largest peak
    
                    scatter(harmonic_ratio, out.rho_max, marker_size, spacing_shapes{s}, 'filled', ...
                            'MarkerFaceColor', colors(w,:), 'markerfacealpha', steepness_alpha(st), ...
                            'HandleVisibility', 'off')
                end
            end
        end
    end

    % Legend
    if c == 1
        % Legend for color
        for w = 1:length(wavelengths)
            plot(nan, nan, 'Color', colors(w,:), 'linewidth', 3, ...
                'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for s = 1:length(farm_spacings)
            scatter(nan, nan, marker_size, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker alpha
        for st = 1:length(wave_steepnesses)
            scatter(nan, nan, marker_size, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
                    'markerfacealpha', steepness_alpha(st), ...
                    'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
        end

        leg = legend('interpreter', 'latex', 'box', 'off');
        leg.Layout.Tile = 'east';
    end
    hold off
end

% Plot axes labels
linkaxes(h, 'xy')
xlim([0.5, 2.5])
ylim([-1, 1])
xlabel(t, '$S_x / \lambda$', 'interpreter', 'latex')
ylabel(t, sprintf('$\\rho_{%s}$', symb), 'interpreter', 'latex')

clear available_wave_cases A B c colors correlating_turbine DOF farm_spacing fixed_row h harmonic_ratio
clear leg maxLagSeconds n name num_waked_turbines out reference_turbine s st steep symb t leg
clear turbine_A_signal turbine_B_signal w waked_turbines wave phase_velocity wave_advection_time


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG: MAX CORRELATION VALUE (WITH SIGN) AGAINST
% HARMOIC RATIO REMAINDER
% LOOPED OVER WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correlating a specified DOF from a fixed row, to all waked rows
% Plotting the absolute value of maximum correlation coefficient

% Which row to correlate against
fixed_row = 1;
DOF = 'pitch_kal';
colors = row_colors.(sprintf('Row%1.0f', fixed_row));


% Cross-correlation
wave_advection_multiplier = 1.5;

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);    

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plotting
clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('Maximum Correlation Coefficient of %s ($%s$): $|\\rho_{%s}|$', name, symb, symb), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    % Get turbines
    reference_turbine = centers(fixed_row);
    correlating_turbine = waked_turbines(c);

    h(c) = nexttile;
    title(sprintf('Row %1.0f to Row %1.0f', fixed_row, ceil(correlating_turbine / 3)))
    hold on
    % Loop through spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through waves
        for st = 1:length(wave_steepnesses)
            steep = compose('%02d', round(100 * wave_steepnesses(st)));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)

                    % How many waves fit between rows
                    harmonic_ratio = farm_spacings(s) / wavelengths(w);
                    harmonic_ratio_remainder = mod(harmonic_ratio, 1);

                    % How long it takes for wave to travel between rows
                    phase_velocity = wavelengths(w) * forcing_frequencies.(['LM', num2str(wavelengths(w))]);
                    wave_advection_time = farm_spacings(s) / phase_velocity;
                    maxLagSeconds = wave_advection_multiplier * wave_advection_time;
    
                    % Get signals
                    turbine_A_signal = tracking.(farm_spacing).(wave)(reference_turbine).(DOF);
                    turbine_B_signal = tracking.(farm_spacing).(wave)(correlating_turbine).(DOF);
                    
                    % Make the same length since we use 'coeff'
                    n = min(numel(turbine_A_signal), numel(turbine_B_signal));
                    A = turbine_A_signal(1:n);
                    B = turbine_B_signal(1:n);
        
                    % Cross-correlate
                    out = xcorr_metrics(A, B, sampling_frequency, maxLagSeconds);
        
                    % out.rho_max ~ largest magnitude XC coefficient
                    % out.rho_abs_max ~ abs of largest magnitude XC coefficient
                    % out.tau_max ~ time lag at largest peak
    
                    scatter(harmonic_ratio_remainder, out.rho_max, marker_size, spacing_shapes{s}, 'filled', ...
                            'MarkerFaceColor', colors(w,:), 'markerfacealpha', steepness_alpha(st), ...
                            'HandleVisibility', 'off')
                end
            end
        end
    end

    % Legend
    if c == 1
        % Legend for color
        for w = 1:length(wavelengths)
            plot(nan, nan, 'Color', colors(w,:), 'linewidth', 3, ...
                'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for s = 1:length(farm_spacings)
            scatter(nan, nan, marker_size, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker alpha
        for st = 1:length(wave_steepnesses)
            scatter(nan, nan, marker_size, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
                    'markerfacealpha', steepness_alpha(st), ...
                    'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
        end

        leg = legend('interpreter', 'latex', 'box', 'off');
        leg.Layout.Tile = 'east';
    end

    % Phase alignment lines
    xline([0, 1], 'linestyle', '--', 'HandleVisibility', 'off', 'label', 'In-Phase', 'fontsize', 8)
    xline(0.5, 'linestyle', '--', 'HandleVisibility', 'off', 'label', 'Out-of-Phase', 'fontsize', 8)
    hold off
end

% Plot axes labels
linkaxes(h, 'xy')
xlim([-0.05, 1.05])
ylim([-1, 1])
xlabel(t, 'Wave Remainder: mod($S_x / \lambda$, 1)', 'interpreter', 'latex')
ylabel(t, sprintf('$\\rho_{%s}$', symb), 'interpreter', 'latex')

clear available_wave_cases A B c colors correlating_turbine DOF farm_spacing fixed_row h harmonic_ratio
clear leg maxLagSeconds n name num_waked_turbines out reference_turbine s st steep symb t leg wave_advection_time
clear turbine_A_signal turbine_B_signal w waked_turbines wave harmonic_ratio_remainder phase_velocity


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG: ABS OF MAX CORRELATION VALUE AGAINST
% HARMOIC RATIO
% LOOPED OVER WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correlating a specified DOF from a fixed row, to all waked rows
% Plotting the absolute value of maximum correlation coefficient

% Which row to correlate against
fixed_row = 1;
DOF = 'yaw_kal';
colors = row_colors.(sprintf('Row%1.0f', fixed_row));


% Cross-correlation
wave_advection_multiplier = 1.5;

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);    

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plotting
clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('ABS of Maximum Correlation Coefficient of %s ($%s$): $|\\rho_{%s}|$', name, symb, symb), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    % Get turbines
    reference_turbine = centers(fixed_row);
    correlating_turbine = waked_turbines(c);

    h(c) = nexttile;
    title(sprintf('Row %1.0f to Row %1.0f', fixed_row, ceil(correlating_turbine / 3)))
    hold on
    % Loop through spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through waves
        for st = 1:length(wave_steepnesses)
            steep = compose('%02d', round(100 * wave_steepnesses(st)));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)

                    % How many waves fit between rows
                    harmonic_ratio = farm_spacings(s) / wavelengths(w);

                    % How long it takes for wave to travel between rows
                    phase_velocity = wavelengths(w) * forcing_frequencies.(['LM', num2str(wavelengths(w))]);
                    wave_advection_time = farm_spacings(s) / phase_velocity;
                    maxLagSeconds = wave_advection_multiplier * wave_advection_time;
    
                    % Get signals
                    turbine_A_signal = tracking.(farm_spacing).(wave)(reference_turbine).(DOF);
                    turbine_B_signal = tracking.(farm_spacing).(wave)(correlating_turbine).(DOF);
                    
                    % Make the same length since we use 'coeff'
                    n = min(numel(turbine_A_signal), numel(turbine_B_signal));
                    A = turbine_A_signal(1:n);
                    B = turbine_B_signal(1:n);
        
                    % Cross-correlate
                    out = xcorr_metrics(A, B, sampling_frequency, maxLagSeconds);
        
                    % out.rho_max ~ largest magnitude XC coefficient
                    % out.rho_abs_max ~ abs of largest magnitude XC coefficient
                    % out.tau_max ~ time lag at largest peak
    
                    scatter(harmonic_ratio, out.rho_abs_max, marker_size, spacing_shapes{s}, 'filled', ...
                            'MarkerFaceColor', colors(w,:), 'markerfacealpha', steepness_alpha(st), ...
                            'HandleVisibility', 'off')
                end
            end
        end
    end

    % Legend
    if c == 1
        % Legend for color
        for w = 1:length(wavelengths)
            plot(nan, nan, 'Color', colors(w,:), 'linewidth', 3, ...
                'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for s = 1:length(farm_spacings)
            scatter(nan, nan, marker_size, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker alpha
        for st = 1:length(wave_steepnesses)
            scatter(nan, nan, marker_size, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
                    'markerfacealpha', steepness_alpha(st), ...
                    'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
        end

        leg = legend('interpreter', 'latex', 'box', 'off');
        leg.Layout.Tile = 'east';
    end
    hold off
end

% Plot axes labels
linkaxes(h, 'xy')
xlim([0.5, 2.5])
ylim([0, 1])
xlabel(t, '$S_x / \lambda$', 'interpreter', 'latex')
ylabel(t, sprintf('$| \\rho_{%s} |$', symb), 'interpreter', 'latex')

clear available_wave_cases A B c colors correlating_turbine DOF farm_spacing fixed_row h harmonic_ratio
clear leg maxLagSeconds n name num_waked_turbines out reference_turbine s st steep symb t leg
clear turbine_A_signal turbine_B_signal w waked_turbines wave wave_advection_time phase_velocity


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG: ABS OF MAX CORRELATION VALUE AGAINST
% HARMOIC RATIO REMAINDER
% LOOPED OVER WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correlating a specified DOF from a fixed row, to all waked rows
% Plotting the absolute value of maximum correlation coefficient

% Which row to correlate against
fixed_row = 1;
DOF = 'yaw_kal';


% Cross-correlation
wave_advection_multiplier = 1.5;

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);    

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plotting
clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('ABS of Maximum Correlation Coefficient of %s ($%s$): $|\\rho_{%s}|$', name, symb, symb), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    % Get turbines
    reference_turbine = centers(fixed_row);
    correlating_turbine = waked_turbines(c);
    colors = row_colors.(sprintf('Row%1.0f', turbineRow(correlating_turbine, farm_arrangement)));

    h(c) = nexttile;
    title(sprintf('Row %1.0f to Row %1.0f', fixed_row, ceil(correlating_turbine / 3)))
    hold on
    % Loop through spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through waves
        for st = 1:length(wave_steepnesses)
            steep = compose('%02d', round(100 * wave_steepnesses(st)));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)

                    % How many waves fit between rows
                    harmonic_ratio = farm_spacings(s) / wavelengths(w);
                    harmonic_ratio_remainder = mod(harmonic_ratio, 1);

                    % How long it takes for wave to travel between rows
                    phase_velocity = wavelengths(w) * forcing_frequencies.(['LM', num2str(wavelengths(w))]);
                    wave_advection_time = farm_spacings(s) / phase_velocity;
                    maxLagSeconds = wave_advection_multiplier * wave_advection_time;
    
                    % Get signals
                    turbine_A_signal = tracking.(farm_spacing).(wave)(reference_turbine).(DOF);
                    turbine_B_signal = tracking.(farm_spacing).(wave)(correlating_turbine).(DOF);
                    
                    % Make the same length since we use 'coeff'
                    n = min(numel(turbine_A_signal), numel(turbine_B_signal));
                    A = turbine_A_signal(1:n);
                    B = turbine_B_signal(1:n);
        
                    % Cross-correlate
                    out = xcorr_metrics(A, B, sampling_frequency, maxLagSeconds);
        
                    % out.rho_max ~ largest magnitude XC coefficient
                    % out.rho_abs_max ~ abs of largest magnitude XC coefficient
                    % out.tau_max ~ time lag at largest peak
    
                    scatter(harmonic_ratio_remainder, out.rho_abs_max, marker_size, spacing_shapes{s}, 'filled', ...
                            'MarkerFaceColor', colors(w,:), 'markerfacealpha', steepness_alpha(st), ...
                            'HandleVisibility', 'off')
                end
            end
        end
    end

    % Legend
    if c == 1
        % Legend for color
        for w = 1:length(wavelengths)
            plot(nan, nan, 'Color', colors(w,:), 'linewidth', 3, ...
                'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for s = 1:length(farm_spacings)
            scatter(nan, nan, marker_size, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker alpha
        for st = 1:length(wave_steepnesses)
            scatter(nan, nan, marker_size, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
                    'markerfacealpha', steepness_alpha(st), ...
                    'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
        end

        leg = legend('interpreter', 'latex', 'box', 'off');
        leg.Layout.Tile = 'east';
    end

    % Phase alignment lines
    xline([0, 1], 'linestyle', '--', 'HandleVisibility', 'off', 'label', 'In-Phase', 'fontsize', 8)
    xline(0.5, 'linestyle', '--', 'HandleVisibility', 'off', 'label', 'Out-of-Phase', 'fontsize', 8)
    hold off
end

% Plot axes labels
linkaxes(h, 'xy')
xlim([-0.05, 1.05])
ylim([0, 1])
xlabel(t, 'Wave Remainder: mod($S_x / \lambda$, 1)', 'interpreter', 'latex')
ylabel(t, sprintf('$| \\rho_{%s} |$', symb), 'interpreter', 'latex')

clear available_wave_cases A B c colors correlating_turbine DOF farm_spacing fixed_row h harmonic_ratio
clear leg maxLagSeconds n name num_waked_turbines out reference_turbine s st steep symb t leg
clear turbine_A_signal turbine_B_signal w waked_turbines wave harmonic_ratio_remainder


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG: TIME LAG AT MAX CORRELATION VALUE AGAINST
% HARMOIC RATIO
% LOOPED OVER WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correlating a specified DOF from a fixed row, to all waked rows
% Plotting the absolute value of maximum correlation coefficient

% Which row to correlate against
fixed_row = 1;
DOF = 'pitch_kal';
colors = row_colors.(sprintf('Row%1.0f', fixed_row));


% Cross-correlation
wave_advection_multiplier = 1.5;

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);    

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plotting
clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('ABS of Maximum Correlation Coefficient of %s ($%s$): $|\\rho_{%s}|$', name, symb, symb), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    % Get turbines
    reference_turbine = centers(fixed_row);
    correlating_turbine = waked_turbines(c);

    h(c) = nexttile;
    title(sprintf('Row %1.0f to Row %1.0f', fixed_row, ceil(correlating_turbine / 3)))
    hold on
    % Loop through spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through waves
        for st = 1:length(wave_steepnesses)
            steep = compose('%02d', round(100 * wave_steepnesses(st)));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)

                    % How many waves fit between rows
                    harmonic_ratio = farm_spacings(s) / wavelengths(w);

                    % How long it takes for wave to travel between rows
                    phase_velocity = wavelengths(w) * forcing_frequencies.(['LM', num2str(wavelengths(w))]);
                    wave_advection_time = farm_spacings(s) / phase_velocity;
                    maxLagSeconds = wave_advection_multiplier * wave_advection_time;
    
                    % Get signals
                    turbine_A_signal = tracking.(farm_spacing).(wave)(reference_turbine).(DOF);
                    turbine_B_signal = tracking.(farm_spacing).(wave)(correlating_turbine).(DOF);
                    
                    % Make the same length since we use 'coeff'
                    n = min(numel(turbine_A_signal), numel(turbine_B_signal));
                    A = turbine_A_signal(1:n);
                    B = turbine_B_signal(1:n);
        
                    % Cross-correlate
                    out = xcorr_metrics(A, B, sampling_frequency, maxLagSeconds);
        
                    % out.rho_max ~ largest magnitude XC coefficient
                    % out.rho_abs_max ~ abs of largest magnitude XC coefficient
                    % out.tau_max ~ time lag at largest peak
    
                    scatter(harmonic_ratio, out.tau_max, marker_size, spacing_shapes{s}, 'filled', ...
                            'MarkerFaceColor', colors(w,:), 'markerfacealpha', steepness_alpha(st), ...
                            'HandleVisibility', 'off')
                end
            end
        end
    end

    % Legend
    if c == 1
        % Legend for color
        for w = 1:length(wavelengths)
            plot(nan, nan, 'Color', colors(w,:), 'linewidth', 3, ...
                'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for s = 1:length(farm_spacings)
            scatter(nan, nan, marker_size, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker alpha
        for st = 1:length(wave_steepnesses)
            scatter(nan, nan, marker_size, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
                    'markerfacealpha', steepness_alpha(st), ...
                    'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
        end

        leg = legend('interpreter', 'latex', 'box', 'off');
        leg.Layout.Tile = 'east';
    end
    hold off
end

% Plot axes labels
linkaxes(h, 'xy')
xlim([0.5, 2.5])
% ylim([0, 1])
xlabel(t, '$S_x / \lambda$', 'interpreter', 'latex')
ylabel(t, sprintf('$\\tau_{\\rho}$ [seconds]'), 'interpreter', 'latex')

clear available_wave_cases A B c colors correlating_turbine DOF farm_spacing fixed_row h harmonic_ratio
clear leg maxLagSeconds n name num_waked_turbines out reference_turbine s st steep symb t leg
clear turbine_A_signal turbine_B_signal w waked_turbines wave




%% Functions

function out = xcorr_metrics(x, y, fs, maxLag_s)
    % Returns peak correlation and lag at peak (seconds)
    % Constrains search to POSITIVE lags only (y lagging x)
    x = x(:); y = y(:);
    x = x - mean(x, 'omitnan');
    y = y - mean(y, 'omitnan');
    
    % Optional: handle NaNs by simple fill (or remove segments)
    x = fillmissing(x,'linear','EndValues','nearest');
    y = fillmissing(y,'linear','EndValues','nearest');
    
    % Perform cross-correlation
    maxLag = round(maxLag_s * fs);
    [r, lags] = xcorr(x, y, maxLag, 'coeff');
    tau = lags / fs;
    
    % Only look at positive lags
    positive_lag_idx = lags >= 0;
    r_positive = r(positive_lag_idx);
    tau_positive = tau(positive_lag_idx);
    
    % % Peak magnitude (use abs so anti-correlation counts as strong coupling)
    % [~, idx] = max(abs(r_positive));
    % out.rho_max = r_positive(idx);
    % out.rho_abs_max = abs(r_positive(idx));
    % out.tau_max = tau_positive(idx);

    % Maximum positive correlation (in-phase coupling)
    [out.rho_max_positive, idx_pos] = max(r_positive);
    out.tau_max_positive = tau_positive(idx_pos);
    
    % Maximum absolute correlation (strongest coupling)
    [out.rho_abs_max, idx_abs] = max(abs(r_positive));
    out.rho_max = r_positive(idx_abs);
    out.tau_max = tau_positive(idx_abs);
    
    % Zero-lag correlation (should always be in positive domain since lag=0 is included)
    [~, iz] = min(abs(tau_positive));
    out.rho0 = r_positive(iz);
end


function row = turbineRow(turbineID, layout)
    % turbineRow: Figures out which row a turbine is in for both inline and
    % staggered configurations

    switch layout
        case 'Inline'
            row = ceil(turbineID / 3);
        case 'Staggered'
            rowSizes = [3 2 3 2];
            rowEdges = [0 cumsum(rowSizes)];
            row = find(turbineID > rowEdges(1:end-1) & ...
                       turbineID <= rowEdges(2:end));
        otherwise
            error('Unknown layout type')
    end
end