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
farm_arrangement = "Inline";
harmonic_cutoff = 5;


% Make fancy title for plots
fancy_title = sprintf('%s Floating Wind Farm', farm_arrangement);

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
freestream_velocity = 4.2;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING PROPERTIES AND CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wavelengths and steepnesses to loop over
wavelengths = [2,3,4,5];
wave_steepnesses = [0.06, 0.09, 0.12];

% Colors per row
row_colors.Row1 = flipud(slanCM(45, 2 * length(farm_spacings)));
row_colors.Row2 = flipud(slanCM(47, 2 * length(farm_spacings)));
row_colors.Row3 = flipud(slanCM(34, 2 * length(farm_spacings)));
row_colors.Row4 = flipud(slanCM(35, 2 * length(farm_spacings)));

% Scatter plot nomenclature
marker_size = 100;
steepness_alpha = [0.3, 0.6, 1];
wave_shapes = {'o', '^', 'square', 'pentagram'};


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ALL CROSS-CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multiplier to either wave/flow advection time to widen cross-correlation 
advection_multiplier = 1.5;
cross_correlations.AdvectionMultiplier = advection_multiplier;

clc; fprintf('Computing cross-correlations...\n')
% Loop through farm arrangements
for a = 1:length(farm_arrangements)
    arrangement = farm_arrangements{a};

    % Get farm layout info
    tmp_centers = turbine_catalog.(farm_arrangement).centers;

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat("WT60_", farm_spacing, "_AG0");
        fprintf('%s: %s\n', arrangement, caze)
        waves = fieldnames(full_tracking.(arrangement).(farm_spacing));
    
        % Loop through waves
        for w = 1:length(waves)

            % Get wave and wave frequency
            wave = waves{w};
            split_wave = split(wave, '_');
            wavelength_name = split_wave{1};

            % Compute wave-advection time between rows [seconds]
            if ~strcmp(wavelength_name, 'LM0')
                wavelength = wavelengths_from_string.(wavelength_name);
                phase_velocity = wavelength * forcing_frequencies.(wavelength_name);
                wave_advection_time_single_row = farm_spacings(s) / phase_velocity;
                % maxLagSeconds = advection_multiplier * wave_advection_time;

                % Save wave advection time to travel across one row of turbines
                wave_advection_times.(farm_spacing).(wave) = wave_advection_time_single_row;

                % Save wave advection time to cross-correlation structure
                cross_correlations.(arrangement).(farm_spacing).(wave).wave_advection_time_single_row = wave_advection_time_single_row;

            % Use a different time scale for no-wave cases
            else
                flow_advection_time = (farm_spacings(s) * 0.15) / freestream_velocity;
                maxLagSeconds = advection_multiplier * flow_advection_time;

                % Save wave advection time to cross-correlation structure
                cross_correlations.(arrangement).(farm_spacing).(wave).flow_advection_time = flow_advection_time;
            end
            

            % Loop through different reference correlating rows (1-3)
            for fixed_row = 1:3

                % Get reference turbine
                reference_turbine = centers(fixed_row);
                reference_row_tag = ['Row', num2str(fixed_row)];

                % Determine waked turbines
                waked_turbines = tmp_centers(fixed_row + 1:end);
                num_waked_turbines = length(waked_turbines);

                % Loop through waked turbines
                for c = 1:num_waked_turbines
                    
                    % Get turbines
                    correlating_turbine = waked_turbines(c);

                    % Max time lag based on time for wave to travel across rows
                    if ~strcmp(wavelength_name, 'LM0')
                        wave_advection_time_across_rows = (c * wave_advection_time_single_row);
                        maxLagSeconds = advection_multiplier * wave_advection_time_across_rows;
                    end

                    % Loop through different DOFs
                    for d = 1:length(DOFs)
                        DOF = DOFs{d};
    
                        % Get signals
                        A = full_tracking.(arrangement).(farm_spacing).(wave)(reference_turbine).(DOF);
                        B = full_tracking.(arrangement).(farm_spacing).(wave)(correlating_turbine).(DOF);
                        
                        % Make the same length since we use 'coeff'
                        n = min(numel(A), numel(B));
                        A = A(1:n);
                        B = B(1:n);
            
                        % Cross-correlate
                        out = xcorr_metrics(A, B, sampling_frequency, maxLagSeconds);

                        % Save
                        cross_correlations.(arrangement).(farm_spacing).(wave).(reference_row_tag)(c).(DOF) = out;
    
                    end
                end
            end
        end
    end
    fprintf('\n')
end

% Save
fprintf('\n')
save_path = fullfile(tracking_path, 'CrossCorrelations', sprintf('OffshoreTracking_PositionCrossCorrelations_SavitskyGolay_Cutoff_%1.0f.mat', harmonic_cutoff));
fprintf('Saving Position Cross-Correlations...\n')
save(save_path, 'cross_correlations')
fprintf('Done Saving!\n\n')

clear a A B c caze correlating_turbine d DOF arrangement farm_spacing fixed_row maxLagSeconds n
clear num_waked_turbines out phase_velocity reference_row_tag reference_turbine s split_wave 
clear tmp_centers w waked_turbines wave wave_advection_time wavelength wavelength_name waves flow_advection_time


%% TESTING ALIGNMENT

spacing = 'SX50';
wave = 'LM5_AK12';

figure('color', 'white')
hold on
plot(tracking.(spacing).(wave)(2).time, tracking.(spacing).(wave)(2).pitch_kal)
plot(tracking.(spacing).(wave)(5).time, tracking.(spacing).(wave)(5).pitch_kal)
hold off
for i = 1:40
    xline(i * (1/forcing_frequencies.LM5), 'linestyle', '--')
end
xlim([0, 20])


%% Hilbert transform test

% Assume x1 and x2 are your two turbine signals, sampled at fs
spacing = 'SX50';
wave = 'LM2_AK12';
f_wave = forcing_frequencies.LM2;
fs = 30;

x1 = tracking.(spacing).(wave)(2).yaw_kal;
x2 = tracking.(spacing).(wave)(5).yaw_kal;
t = tracking.(spacing).(wave)(2).time;

% 1. Bandpass around your dominant frequency to clean things up
bw = 0.5;    % bandwidth (adjust as needed)
x1_filt = bandpass(x1, [f_wave - bw, f_wave + bw], fs);
x2_filt = bandpass(x2, [f_wave - bw, f_wave + bw], fs);
% x1_filt = x1;
% x2_filt = x2;

% 2. Compute analytic signals
z1 = hilbert(x1_filt);
z2 = hilbert(x2_filt);

% 3. Extract instantaneous phase
phi1 = unwrap(angle(z1));
phi2 = unwrap(angle(z2));

% 4. Phase difference
delta_phi = phi2 - phi1;

% 5. Convert to time lag at the dominant frequency
% Phase in radians, frequency in Hz -> lag in seconds
lag_instantaneous = delta_phi / (2 * pi * f_wave);

% 6. Visualize
% t = (0:length(x1)-1) / fs;

figure
subplot(3,1,1)
plot(t, x1_filt, t, x2_filt)
xlabel('Time (s)')
ylabel('Filtered signal')
legend('Turbine 1', 'Turbine 2')

subplot(3,1,2)
plot(t, delta_phi)
xlabel('Time (s)')
ylabel('\Delta\phi (rad)')
title('Instantaneous phase difference')

subplot(3,1,3)
plot(t, lag_instantaneous)
xlabel('Time (s)')
ylabel('Equivalent lag (s)')
title('Instantaneous time lag')

% Mean phase difference (in radians)
mean_delta_phi = mean(delta_phi);

% Convert to time lag
mean_lag = mean_delta_phi / (2 * pi * f_wave);

% You can also look at the distribution
figure
histogram(delta_phi, 50)
xlabel('\Delta\phi (rad)')
ylabel('Count')
title('Distribution of instantaneous phase difference')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG: MAX CORRELATION VALUE (WITH SIGN) AGAINST
% WAVELENGTH
% LOOPED OVER WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rho_pos_max ~ max positive correlation (positive definite)
% tau_pos_max ~ time lag associated with max positive correlation
% rho_abs_max ~ max magnitude correlation (negative or positive)
% tau_abs_max ~ time lag associated with max magnitude correlation

% Which row to correlate against
fixed_row = 1;
DOF = 'pitch_kal';
cross_correlation_metric = 'rho_pos_max';


% Name fixed row
reference_row_tag = ['Row', num2str(fixed_row)];
colors = row_colors.(reference_row_tag);

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF); 

% Adjust plot based on value plotted
plotDetails = crossCorrelationPlotDetails(cross_correlation_metric, symb);

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plotting
clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('%s\n%s of %s ($%s$): %s', fancy_title, plotDetails.title, name, symb, plotDetails.symb), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    % Get turbines
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

                    % Load cross-correlation metric
                    data = cross_correlations.(farm_arrangement).(farm_spacing).(wave).(reference_row_tag)(c).(DOF).(cross_correlation_metric);
        
                    % Plot
                    scatter(wavelengths(w), data, marker_size, 'filled', ...
                            'MarkerFaceColor', colors(s,:), ...
                            'MarkerFaceAlpha', steepness_alpha(st), ...
                            'Marker', wave_shapes{w}, ...
                            'HandleVisibility', 'off')

                    tmpX(w) = wavelengths(w);
                    tmpY(w) = data;
                end
            end
            % Connect the dots
            P = plot(tmpX, tmpY, 'color', colors(s,:), 'linewidth', 1, 'HandleVisibility', 'off');
            P.Color(4) = steepness_alpha(st);
        end
    end

    % Legend
    if c == 1
        % Legend for color
        for s = 1:length(farm_spacings)
            plot(nan, nan, 'Color', colors(s,:), 'linewidth', 3, ...
                'Displayname', sprintf('$S_x = %1.1fD', farm_spacings(s)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for w = 1:length(wavelengths)
            scatter(nan, nan, marker_size, wave_shapes{w}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$\\lambda = %1.0fD$', wavelengths(w)))
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
    xticks(2:1:5)
end

% Plot axes labels
linkaxes(h, 'xy')
xlim([1.5, 5.5])
ylim(plotDetails.ylims)
xlabel(t, '$\lambda / D$', 'interpreter', 'latex')
ylabel(t, plotDetails.symb, 'interpreter', 'latex')

clear available_wave_cases A B c colors correlating_turbine DOF farm_spacing fixed_row h harmonic_ratio
clear leg maxLagSeconds n name num_waked_turbines out reference_turbine s st steep symb t leg reference_row_tag plotDetails
clear turbine_A_signal turbine_B_signal w waked_turbines wave phase_velocity wave_advection_time data cross_correlation_metric


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG: MAX CORRELATION VALUE (WITH SIGN) AGAINST
% HARMOIC RATIO
% LOOPED OVER WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rho_pos_max ~ max positive correlation (positive definite)
% tau_pos_max ~ time lag associated with max positive correlation
% rho_abs_max ~ max magnitude correlation (negative or positive)
% tau_abs_max ~ time lag associated with max magnitude correlation

% Which row to correlate against
fixed_row = 1;
DOF = 'pitch_kal';
cross_correlation_metric = 'rho_pos_max';


% Name fixed row
reference_row_tag = ['Row', num2str(fixed_row)];
colors = row_colors.(reference_row_tag);

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF); 

% Adjust plot based on value plotted
plotDetails = crossCorrelationPlotDetails(cross_correlation_metric, symb);

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plotting
clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('%s\n%s of %s ($%s$): %s', fancy_title, plotDetails.title, name, symb, plotDetails.symb), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    % Get turbines
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

                    % Load cross-correlation metric
                    data = cross_correlations.(farm_arrangement).(farm_spacing).(wave).(reference_row_tag)(c).(DOF).(cross_correlation_metric);
        
                    % Plot
                    scatter(harmonic_ratio, data, marker_size, 'filled', ...
                            'MarkerFaceColor', colors(s,:), ...
                            'MarkerFaceAlpha', steepness_alpha(st), ...
                            'Marker', wave_shapes{w}, ...
                            'HandleVisibility', 'off')
                end
            end
        end
    end

    % Legend
    if c == 1
        % Legend for color
        for s = 1:length(farm_spacings)
            plot(nan, nan, 'Color', colors(s,:), 'linewidth', 3, ...
                'Displayname', sprintf('$S_x = %1.1fD', farm_spacings(s)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for w = 1:length(wavelengths)
            scatter(nan, nan, marker_size, wave_shapes{w}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$\\lambda = %1.0fD$', wavelengths(w)))
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
ylim(plotDetails.ylims)
xlabel(t, '$S_x / \lambda$', 'interpreter', 'latex')
ylabel(t, plotDetails.symb, 'interpreter', 'latex')

clear available_wave_cases A B c colors correlating_turbine DOF farm_spacing fixed_row h harmonic_ratio
clear leg maxLagSeconds n name num_waked_turbines out reference_turbine s st steep symb t leg reference_row_tag plotDetails
clear turbine_A_signal turbine_B_signal w waked_turbines wave phase_velocity wave_advection_time data cross_correlation_metric


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG: MAX CORRELATION VALUE TIME LAG (ALIGNEMNT TIME) 
% AGAINST WAVE PHASE MISALIGNMENT
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rho_pos_max ~ max positive correlation (positive definite)
% tau_pos_max ~ time lag associated with max positive correlation
% rho_abs_max ~ max magnitude correlation (negative or positive)
% tau_abs_max ~ time lag associated with max magnitude correlation

% Which row to correlate against
fixed_row = 1;
DOF = 'pitch_kal';
cross_correlation_metric = 'tau_pos_max';


% Name fixed row
reference_row_tag = ['Row', num2str(fixed_row)];
colors = row_colors.(reference_row_tag);

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF); 

% Adjust plot based on value plotted
plotDetails = crossCorrelationPlotDetails(cross_correlation_metric, symb);

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plotting
clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('%s\n%s of %s ($%s$): %s', fancy_title, plotDetails.title, name, symb, plotDetails.symb), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    % Get turbines
    correlating_turbine = waked_turbines(c);

    h(c) = nexttile;
    title(sprintf('Row %1.0f to Row %1.0f', fixed_row, ceil(correlating_turbine / 3)))
    hold on
    % Loop through spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through waves
        % for st = 1:length(wave_steepnesses)
        for st = 3
            steep = compose('%02d', round(100 * wave_steepnesses(st)));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Time for wave to travel from one row to another
                wave_advection_time = wave_advection_times.(farm_spacing).(wave);

                % Wave period time scale
                wave_period = 1 / forcing_frequencies.(['LM', num2str(wavelengths(w))]);

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)

                    % Load cross-correlation metric
                    data = cross_correlations.(farm_arrangement).(farm_spacing).(wave).(reference_row_tag)(c).(DOF).(cross_correlation_metric);

                    % How many waves fit between rows
                    distance_between_rows = c * farm_spacings(s);
                    wavelengths_between_turbines = distance_between_rows / wavelengths(w);
                    wave_phase_misalignment = mod(wavelengths_between_turbines, 1);

                    
                    % Wave travel time between rows time scale
                    wave_advection_time_between_rows = c * wave_advection_time;

                    % Wake advection time
                    wake_advection_time_between_rows = (c * farm_spacings(s) * 0.15) / freestream_velocity;

                    data = mod(data, 1);

        
                    % Plot
                    scatter(farm_spacings(s), data, marker_size, 'filled', ...
                            'MarkerFaceColor', colors(s,:), ...
                            'MarkerFaceAlpha', steepness_alpha(st), ...
                            'Marker', wave_shapes{w}, ...
                            'HandleVisibility', 'off')
                end
            end
        end
    end

    % Legend
    if c == 1
        % Legend for color
        for s = 1:length(farm_spacings)
            plot(nan, nan, 'Color', colors(s,:), 'linewidth', 3, ...
                'Displayname', sprintf('$S_x = %1.1fD', farm_spacings(s)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for w = 1:length(wavelengths)
            scatter(nan, nan, marker_size, wave_shapes{w}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$\\lambda = %1.0fD$', wavelengths(w)))
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
    xticks(2:1:5)
end

% Plot axes labels
linkaxes(h, 'xy')
xlim([1.5, 5.5])
% ylim(plotDetails.ylims)
xlabel(t, '$\lambda / D$', 'interpreter', 'latex')
ylabel(t, plotDetails.symb, 'interpreter', 'latex')

clear available_wave_cases A B c colors correlating_turbine DOF farm_spacing fixed_row h harmonic_ratio
clear leg maxLagSeconds n name num_waked_turbines out reference_turbine s st steep symb t leg reference_row_tag plotDetails
clear turbine_A_signal turbine_B_signal w waked_turbines wave phase_velocity wave_advection_time data cross_correlation_metric



%% Functions

function out = crossCorrelationPlotDetails(metric, DOF_symb)

    % Correlation coefficient latex symbols
    if contains(metric, 'rho_pos') 
        out.symb = sprintf('$\\max{ \\left( \\rho_{%s} \\right)}$', DOF_symb);
        out.title = 'Maximum Positive Cross-Correlation Value';
        out.ylims = [0,1];
    end
    if contains(metric, 'rho_abs') 
        out.symb = sprintf('$\\max{ \\left( |\\rho_{%s}| \\right)}$', DOF_symb);
        out.title = 'Maximum Magnitude Cross-Correlation Value';
        out.ylims = [-1,1];
    end

    % Tau latex symbols
    if contains(metric, 'tau') && contains(metric, 'pos') 
        out.symb = sprintf('$\\tau_{%s}$', DOF_symb);
        out.title = 'Time Lag at Maximum Positive Cross-Correlation Value';
    end
    if contains(metric, 'tau') && contains(metric, 'abs') 
        out.symb = sprintf('$\\tau_{%s}$', DOF_symb);
        out.title = 'Time Lag at Maximum Magnitude Cross-Correlation Value';
    end

end


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
    r = r(positive_lag_idx);
    tau = tau(positive_lag_idx);
    
    % Save cross-correlation signals and max lag
    out.rho_signal = r;
    out.tau_signal = tau;
    out.input_max_lag_seconds = maxLag_s;

    % Maximum positive correlation (in-phase coupling)
    [~, idx_pos] = max(r);
    out.rho_pos_max = r(idx_pos);
    out.tau_pos_max = tau(idx_pos);
    
    % Maximum absolute correlation (strongest coupling)
    [~, idx_abs] = max(abs(r));
    out.rho_abs_max = r(idx_abs);
    out.tau_abs_max = tau(idx_abs);
    
    % Zero-lag correlation (should always be in positive domain since lag=0 is included)
    [~, iz] = min(abs(tau));
    out.rho_no_lag = r(iz);
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