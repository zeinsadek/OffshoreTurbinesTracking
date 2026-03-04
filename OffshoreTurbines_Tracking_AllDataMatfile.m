%% Combining all the tracking data into a single matfile
% Also filtering and removing obvious outliers and identifying problematic
% cases

% Zein Sadek
% 1/8/2025

clear; close all; clc;
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT TRACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

offshore_path = "/Users/zeinsadek/Desktop/Experiments/Offshore";
tracking_path = fullfile(offshore_path, "Tracking/Data/Matfiles");

% Load both farm arrangements
farm_arrangements = {'Inline', 'Staggered'};

% Load all spacings [D]
farm_spacings = [5, 4.5, 4, 3.5, 3];

% Loop through all farm layouts
clc;
for a = 1:length(farm_arrangements)
    farm_arrangement = farm_arrangements{a};
    fprintf('Loading %s Tracking Data...\n\n', farm_arrangement)

    % Look through all different farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat("WT60_", farm_spacing, "_AG0");
        fprintf('%s: %s\n', farm_arrangement, caze)
        
        % Import Tracking
        tracking_file = fullfile(tracking_path, farm_arrangement, strcat(caze, ".mat"));
        tmp = load(tracking_file);
        tmp = tmp.(caze);
        tmp_tracking.(farm_arrangement).(farm_spacing) = tmp;
    end
    fprintf('\n')
end

clear s a caze tracking_file tmp farm_spacing farm_arrangement



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODIFY STRUCTURE FIELD NAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert fieldnames into arry for turbines in tracking structure to match
% power
% Index turbines instead of having them as fields

% Loop through farm layouts
clc;
for a = 1:length(farm_arrangements)
    farm_arrangement = farm_arrangements{a};
    fprintf('Restructuring %s Turbine Field Names...\n\n', farm_arrangement)

    % Loop through different farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat(farm_arrangement, ": WT60_", farm_spacing, "_AG0");
        fprintf('%s\n', caze)
        tmp = tmp_tracking.(farm_arrangement).(farm_spacing);
    
        waves = fieldnames(tmp);
        turbines = fieldnames(tmp.(waves{1}));

        % Save turbine names for both arrangements
        turbine_catalog.(farm_arrangement) = turbines;
    
        % Loop through different waves
        for w = 1:length(waves)
            wave = waves{w};

            % Loop through different turbines
            for t = 1:length(turbines)
                turbine = turbines{t};
                tracking.(farm_arrangement).(farm_spacing).(wave)(t) = tmp_tracking.(farm_arrangement).(farm_spacing).(wave).(turbine);
            end
        end
    end
    fprintf('\n')
end

clear s w t a wave turbine
clear farm_spacing caze tmp waves



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECT COORDINATE SYSTEM + FILL NANS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correct coordinate system to be aligned with:
% +x - streamwise direction
% +y - wall-normal direction
% +z - pointing towards wind-tunnel doors

% Fill in any NaNs that may be in the signals
% Convert translations in meters from centimeters
% Rotations are saved in degrees

fillmethod = 'spline';

% Loop through different farm arrangements
clc;
for a = 1:length(farm_arrangements)
    farm_arrangement = farm_arrangements{a};
    fprintf('Correcting %s Farm Coordinate System...\n\n', farm_arrangement)

    % Loop through different spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat(farm_arrangement, ": WT60_", farm_spacing, "_AG0");
        fprintf('%s\n', caze)
        
        tmp = tracking.(farm_arrangement).(farm_spacing);
        waves = fieldnames(tmp);
    
        % Loop through different waves
        for w = 1:length(waves)
            wave = waves{w};

            % Loop through different turbines
            for t = 1:length(turbine_catalog.(farm_arrangement))
        
                % Time
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).time = tracking.(farm_arrangement).(farm_spacing).(wave)(t).time;
        
                %%% Raw
                % Convert cm to m
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).x = fillmissing(tracking.(farm_arrangement).(farm_spacing).(wave)(t).x, fillmethod) .* 1E-2;
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).y = fillmissing(tracking.(farm_arrangement).(farm_spacing).(wave)(t).z, fillmethod) .* 1E-2;
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).z = fillmissing(-1 * tracking.(farm_arrangement).(farm_spacing).(wave)(t).y, fillmethod) .* 1E-2;
        
                % Rotation in degrees
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).roll = fillmissing(tracking.(farm_arrangement).(farm_spacing).(wave)(t).roll, fillmethod);
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).pitch = fillmissing(-1 * tracking.(farm_arrangement).(farm_spacing).(wave)(t).pitch, fillmethod);
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).yaw = fillmissing(tracking.(farm_arrangement).(farm_spacing).(wave)(t).yaw, fillmethod);
                
                %%% Kalman
                % Convert cm to m
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).x_kal = fillmissing(tracking.(farm_arrangement).(farm_spacing).(wave)(t).x_kal, fillmethod) .* 1E-2;
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).y_kal = fillmissing(tracking.(farm_arrangement).(farm_spacing).(wave)(t).z_kal, fillmethod) .* 1E-2;
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).z_kal = fillmissing(-1 * tracking.(farm_arrangement).(farm_spacing).(wave)(t).y_kal, fillmethod) .* 1E-2;
        
                % Rotation in degrees
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).roll_kal = fillmissing(tracking.(farm_arrangement).(farm_spacing).(wave)(t).roll_kal, fillmethod);
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).pitch_kal = fillmissing(-1 * tracking.(farm_arrangement).(farm_spacing).(wave)(t).pitch_kal, fillmethod);
                corrected_tracking.(farm_arrangement).(farm_spacing).(wave)(t).yaw_kal = fillmissing(tracking.(farm_arrangement).(farm_spacing).(wave)(t).yaw_kal, fillmethod);
        
            end
        end
    end
    fprintf('\n')
end

tracking = corrected_tracking;


% Save an unfiltered version for trackable errors
save_path = fullfile(tracking_path, "AllData", "OffshoreTracking_AllDataCombined_NoFilters_AllCases.mat");
clc; fprintf('Saving Combined, No Filter Data...\n')
save(save_path, 'tracking')
fprintf('Done Saving!\n\n')

clear s w t a caze corrected_tracking farm_spacing save_path
clear fillmethod tmp tmp_tracking wave waves turbines farm_arrangement



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTERING DATA AND IDENTIFYING PROBLEMATIC CASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Look at signals and identify obvious outliers and non-physical
% measurements

% Degrees of freedom
DOFs = {'x_kal', 'y_kal', 'z_kal', 'roll_kal', 'pitch_kal', 'yaw_kal'};

% Limit to define a case as problematic
std_lim = 100;

% Counters
signal_counter = 1;
problems = 1;

% Loop through different farm arrangements
clc;
for a = 1:length(farm_arrangements)
    farm_arrangement = farm_arrangements{a};
    fprintf('Inspecting %s Farm Signals...\n\n', farm_arrangement)

    % Loop through different spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat(farm_arrangement, ": WT60_", farm_spacing, "_AG0");
        fprintf('%s\n', caze)
        
        tmp = tracking.(farm_arrangement).(farm_spacing);
        waves = fieldnames(tmp);
            
        % Loop through different waves
        for w = 1:length(waves)
            wave = waves{w};

            % Loop through different turbines
            for t = 1:length(turbine_catalog.(farm_arrangement))
                
                % Loop through degrees-of-freedom
                for d = 1:length(DOFs)
                    DOF = DOFs{d};

                    % Compute STD of signals and flag if too large
                    stddev = std(tracking.(farm_arrangement).(farm_spacing).(wave)(t).(DOF), 0, 'all', 'omitnan');
                    devs(signal_counter) = stddev;

                    if stddev > std_lim
                        fprintf('%s, Turbine %1.0f, %s %s\n', caze, t, wave, DOF)
                        problematic_cases(problems).arrangement = farm_arrangement;
                        problematic_cases(problems).spacing = farm_spacing;
                        problematic_cases(problems).wave = wave;
                        problematic_cases(problems).turbine = t;
                        problematic_cases(problems).DOF = DOF;

                        problems = problems + 1;
                    end
                    signal_counter = signal_counter + 1;
                end
            end
        end
    end
end

% Fix extra count at end of loop
problems = problems - 1;
fprintf('\n%2.0f Problematic cases out of %4.0f total cases (%2.2f percent)\n\n', problems, signal_counter, (problems / signal_counter) * 100)

figure('color', 'white')
plot(devs)
title('Standard Deviation of All Tracking Signals')
ylabel('Standard Deviation')
xlabel('Signal')
yline(std_lim, 'linestyle', '--')
yscale('log')

clear a signal_counter caze d DOF farm_arrangement farm_spacing
clear s t tmp w wave waves stddev std_lim devs



%% Plot the problematic cases

% clc; close all
% for p = 1:problems
%     farm_arrangement = problemtatic_cases(p).arrangement;
%     farm_spacing = problemtatic_cases(p).spacing;
%     wave = problemtatic_cases(p).wave;
%     turbine = problemtatic_cases(p).turbine;
%     DOF = problemtatic_cases(p).DOF;
% 
%     caze = strcat(farm_arrangement, ": WT60_", farm_spacing, "_AG0");
%     figure('color', 'white')
%     hold of
%     plot(tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).time, tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).(DOF))
%     plot(tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).time, tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).(DOF(1:end - 4)))
%     hold off
%     title(sprintf('%s, Turbine %1.0f, %s', caze, turbine, DOF), 'interpreter', 'none')
% end
% 
% clear p farm_arrangement farm_spacing wave turbine DOF caze 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE PROBLEMATIC CASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specifically we are removing:
% Staggered, Sx = 3D, LM5_AK12 case entirely
% Too many turbines had entirely incorrect readings to trust the data set

tracking_fields = {'time', 'x', 'y', 'z', 'roll' ,'pitch', 'yaw', ...
                   'x_kal', 'y_kal', 'z_kal', 'roll_kal', 'pitch_kal', 'yaw_kal'};

% Loop through different farm arrangements
clc;
for a = 1:length(farm_arrangements)
    farm_arrangement = farm_arrangements{a};
    fprintf('Correcting %s Farm Coordinate System...\n\n', farm_arrangement)

    % Loop through different spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat(farm_arrangement, ": WT60_", farm_spacing, "_AG0");
        fprintf('%s\n', caze)
        
        tmp = tracking.(farm_arrangement).(farm_spacing);
        waves = fieldnames(tmp);
    
        % Loop through different waves
        for w = 1:length(waves)
            wave = waves{w};

            % Check for problematic case here
            BadCase = strcmp(farm_arrangement, 'Staggered') && ...
                      strcmp(farm_spacing, 'SX30') && ...
                      strcmp(wave, 'LM5_AK12');
            if BadCase
                continue
            else
                % Loop through different turbines
                for t = 1:length(turbine_catalog.(farm_arrangement))
                    
                    % Loop through tracking fields
                    for d = 1:length(tracking_fields)
                        tracking_field = tracking_fields{d};
                        non_problematic_tracking.(farm_arrangement).(farm_spacing).(wave)(t).(tracking_field) = tracking.(farm_arrangement).(farm_spacing).(wave)(t).(tracking_field);
                    end

                end
            end
        end
    end
    fprintf('\n')
end

clear s t tmp tracking_field turbine w wave waves
clear a caze d farm_arrangement farm_spacing p problems
clear problematic_cases c BadCase tracking_fields problematic_cases



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK STATISTICS AND SAVE VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Degrees of freedom
DOFs = {'x_kal', 'y_kal', 'z_kal', 'roll_kal', 'pitch_kal', 'yaw_kal'};

% Clear problem trackers
clear problems signal_counter problematic_cases

% Limit to define a case as problematic
std_lim = 100;

% Counters
signal_counter = 1;
problems = 1;

% Loop through different farm arrangements
clc;
for a = 1:length(farm_arrangements)
    farm_arrangement = farm_arrangements{a};
    fprintf('Inspecting %s Farm Signals...\n\n', farm_arrangement)

    % Loop through different spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat(farm_arrangement, ": WT60_", farm_spacing, "_AG0");
        fprintf('%s\n', caze)
        
        tmp = non_problematic_tracking.(farm_arrangement).(farm_spacing);
        waves = fieldnames(tmp);
            
        % Loop through different waves
        for w = 1:length(waves)
            wave = waves{w};

            % Loop through different turbines
            for t = 1:length(turbine_catalog.(farm_arrangement))
                
                % Loop through degrees-of-freedom
                for d = 1:length(DOFs)
                    DOF = DOFs{d};

                    % Compute STD of signals and flag if too large
                    signal = non_problematic_tracking.(farm_arrangement).(farm_spacing).(wave)(t).(DOF);
                    stddev = std(signal, 0, 'all', 'omitnan');
                    devs(signal_counter) = stddev;

                    % Check for nans
                    if sum(isnan(signal)) > 0
                        fprintf('%s, Turbine %1.0f, %s %s has NaNs\n', caze, t, wave, DOF)
                    end

                    if stddev > std_lim
                        fprintf('%s, Turbine %1.0f, %s %s\n', caze, t, wave, DOF)
                        problematic_cases(problems).arrangement = farm_arrangement;
                        problematic_cases(problems).spacing = farm_spacing;
                        problematic_cases(problems).wave = wave;
                        problematic_cases(problems).turbine = t;
                        problematic_cases(problems).DOF = DOF;

                        problems = problems + 1;
                    end
                    signal_counter = signal_counter + 1;
                end
            end
        end
    end
end

% Rename the main working variable
tracking = non_problematic_tracking;

% Fix extra count at end of loop
problems = problems - 1;
fprintf('\n%2.0f Problematic cases out of %4.0f total cases (%2.2f percent)\n\n', problems, signal_counter, (problems / signal_counter) * 100)

figure('color', 'white')
tiledlayout(2,1)
nexttile
plot(devs)
title('Standard Deviation of All Tracking Signals')
ylabel('Standard Deviation')
xlabel('Signal')
yline(std_lim, 'linestyle', '--')
xlim([0, signal_counter])
yscale('log')

nexttile
hist(devs,20)
fprintf('Average STD: %2.2f\nMax STD %2.2f\n\n', mean(devs, 'all', 'omitnan'), max(devs, [], 'all', 'omitnan'))

% Save a version with the problematic cases removed
save_path = fullfile(tracking_path, "AllData", "OffshoreTracking_AllDataCombined_NoFilters_ProblematicCasesRemoved.mat");
fprintf('Saving Combined, No Filter, Non-Problematic Data...\n')
save(save_path, 'tracking')
fprintf('Done Saving!\n\n')

clear a signal_counter caze d DOF farm_arrangement farm_spacing problems devs save_path
clear s t tmp w wave waves stddev std_lim non_problematic_tracking problematic cases signal



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVITSKY-GOLAY FILTER REMAINING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Testing different cutoff frequencies to apply to all signals, as a
% function of wave forcing frequency

% Sampling frequency [Hz]
sampling_frequency = 30;

% Forced wave frequencies [Hz]
forcing_frequencies.("LM5")  = 1.4273;
forcing_frequencies.("LM4")  = 1.6075;
forcing_frequencies.("LM33") = 1.7651;
forcing_frequencies.("LM3")  = 1.8617;
forcing_frequencies.("LM25") = 2.0402;
forcing_frequencies.("LM2")  = 2.2813;

% Test case
farm_arrangement = 'Inline';
farm_spacing = 'SX50';
wave = 'LM5_AK12';
turbine = 2;

% Get wavelength and frequency
tmp = split(wave, '_');
wavelength = tmp{1};
wave_frequency = forcing_frequencies.(wavelength);
clear tmp

% Get sgolay parameters
harmonic_cutoff = 4;
[polynomial_order, framelength] = sgolay_params(wave_frequency, sampling_frequency, harmonic_cutoff);
clc; fprintf('Harmonic cutoff = %1.1f\nPolynomial Order = %1.0f\nWindow length = %2.0f\n\n', harmonic_cutoff, polynomial_order, framelength)


% Plot to compare
close all
figure('color', 'white')
tiledlayout(2,3)
sgtitle(sprintf('Cutoff Frequency = $%1.0f \\cdot f_{wave}$', harmonic_cutoff), 'interpreter', 'latex')
for d = 1:length(DOFs)
    DOF = DOFs{d};
    h(d) = nexttile;
    hold on
    % Original signal
    plot(tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).time, tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).(DOF), ...
         'linewidth', 2, 'Color', 'red', 'displayname', 'Original')

    % Sgolay-filt signal
    plot(tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).time, sgolayfilt(tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).(DOF), polynomial_order, framelength), ...
         'linewidth', 2, 'Color', 'blue', 'displayname', 'Repaired')
    hold off
    xlim([0, 5])
    legend()
    title(DOF, 'interpreter', 'none')
end

linkaxes(h, 'x')

clear h harmonic_cutoff steepness framelength polynomial_order turbine wave_frequency d DOF
clear farm_spacing wave wavelength farm_arrangement 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVITSKY-GOLAY FILTER REMAINING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Applying a Savitsky-Golay filter to all the Kalman signals
% Polynomial order is set to 3 and filter length is a function of wave
% frequency
% Cuttoff frequency is varied to get different levels of smoothing


% Cutoff frequency multiple of wave frequency
harmonic_cutoff = 5;

% Clears space
clear devs

% Counters
signal_counter = 1;

% Loop through different farm arrangements
clc;
for a = 1:length(farm_arrangements)
    farm_arrangement = farm_arrangements{a};
    fprintf('Filtering %s Farm...\n\n', farm_arrangement)

    % Loop through different spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat(farm_arrangement, ": WT60_", farm_spacing, "_AG0");
        fprintf('%s\n', caze)
        
        tmp = tracking.(farm_arrangement).(farm_spacing);
        waves = fieldnames(tmp);
    
        % Loop through different waves
        for w = 1:length(waves)
            wave = waves{w};

            % Loop through different turbines
            for t = 1:length(turbine_catalog.(farm_arrangement))

                % Add time vector to structure
                filtered_tracking.(farm_arrangement).(farm_spacing).(wave)(t).time = tracking.(farm_arrangement).(farm_spacing).(wave)(t).time;
                
                % Loop through tracking fields
                for d = 1:length(DOFs)
                    DOF = DOFs{d};

                    % Get signal
                    signal = tracking.(farm_arrangement).(farm_spacing).(wave)(t).(DOF);
                    
                    % Get wavelength and frequency
                    tmp = split(wave, '_');
                    wavelength = tmp{1};
                    clear tmp

                    % If wave case, apply Savitsky-Golay filter
                    if ~strcmp(wavelength, 'LM0')
                        % Get wave frequency
                        wave_frequency = forcing_frequencies.(wavelength);
    
                        % Get sgolay parameters
                        [polynomial_order, framelength] = sgolay_params(wave_frequency, sampling_frequency, harmonic_cutoff);
    
                        % Filter signal
                        filtered_signal = sgolayfilt(tracking.(farm_arrangement).(farm_spacing).(wave)(t).(DOF), polynomial_order, framelength);
                        
                        % Save repaired signal
                        filtered_tracking.(farm_arrangement).(farm_spacing).(wave)(t).(DOF) = filtered_signal;

                        % Save filter design to structure
                        filtered_tracking.FilterDesign.(wave).order = polynomial_order;
                        filtered_tracking.FilterDesign.(wave).length = framelength;
                        filtered_tracking.FilterDesign.(wave).cutoff = harmonic_cutoff;

                    % If no-wave case, skip filtering
                    else
                        filtered_tracking.(farm_arrangement).(farm_spacing).(wave)(t).(DOF) = signal;
                    end

                    % Compute STD of signals and flag if too large
                    stddev = std(filtered_tracking.(farm_arrangement).(farm_spacing).(wave)(t).(DOF), 0, 'all', 'omitnan');
                    devs(signal_counter) = stddev;

                    signal_counter = signal_counter + 1;
                end

            end
        end
    end
    fprintf('\n')
end



figure('color', 'white')
tiledlayout(2,1)
nexttile
plot(devs)
title('Standard Deviation of All Tracking Signals')
ylabel('Standard Deviation')
xlabel('Signal')
xlim([0, signal_counter])
yscale('log')

nexttile
hist(devs,20)
clc; fprintf('Average STD: %2.2f\nMax STD %2.2f\n\n', mean(devs, 'all', 'omitnan'), max(devs, [], 'all', 'omitnan'))

% Rename main working variable
tracking = filtered_tracking;

% Save a filtered kalman signals with cutoff frequency in filename
save_path = fullfile(tracking_path, "AllData", sprintf("OffshoreTracking_AllDataCombined_SavitskyGolay_Cutoff_%1.0f.mat", harmonic_cutoff));
fprintf('Saving Combined, Savitsky-Golay Filtered Data...\n')
save(save_path, 'tracking')
fprintf('Done Saving!\n\n')

clear a caze d devs DOF farm_arrangement farm_spacing filtered_signal framelength polynomial_order
clear s signal signal_counter stddev t w wave wave_frequency wavelength waves







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [order, framelen] = sgolay_params(f_drive, fs, cutoff_multiplier)
% SGOLAY_PARAMS Select Savitzky-Golay parameters based on driving frequency
%
%   [order, framelen] = sgolay_params(f_drive, fs)
%   [order, framelen] = sgolay_params(f_drive, fs, cutoff_multiplier)
%
%   Inputs:
%       f_drive           - Driving/wave frequency (Hz)
%       fs                - Sampling frequency (Hz)
%       cutoff_multiplier - Cutoff = multiplier * f_drive (default: 4)
%
%   Outputs:
%       order    - Polynomial order
%       framelen - Frame length (odd integer)

    if nargin < 3, cutoff_multiplier = 4; end
    
    % Target cutoff frequency
    f_cutoff = cutoff_multiplier * f_drive;
    
    % Polynomial order (3 is a good default for smoothing while preserving shape)
    order = 3;
    
    % Compute frame length from cutoff frequency
    framelen = round((order + 1) * fs / (3.2 * f_cutoff));
    
    % Ensure odd
    if mod(framelen, 2) == 0
        framelen = framelen + 1;
    end
    
    % Ensure framelen > order
    if framelen <= order
        framelen = order + 2;
        if mod(framelen, 2) == 0
            framelen = framelen + 1;
        end
    end
end






