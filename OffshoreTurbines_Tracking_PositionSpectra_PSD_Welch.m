%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFFSHORE TURBINES: PSD OF DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at power spectral density of motion across floating wind farm
% Code plots PSD of different degrees-of-freedom and saves matfile with all
% values, with units and normalized by wave-size

% Zein Sadek
% 1/25

clear; close all; clc;
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT TRACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Paths
offshore_path = "/Users/zeinsadek/Desktop/Experiments/Offshore";
tracking_path = fullfile(offshore_path, "Tracking/Data/Matfiles");

% Farm layout + Data filter version
farm_arrangement = "Inline";
harmonic_cutoff = 5;


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
% Farm spacings
spacings_from_string.("SX50") = 5;
spacings_from_string.("SX45") = 4.5;
spacings_from_string.("SX40") = 4;
spacings_from_string.("SX35") = 3.5;
spacings_from_string.("SX30") = 3;
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
% PSD COMPUTATION OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% \Delta f = f_{wave} / n_periods

% n_period controls the frequency resolution and also smoothess of PSD
% Smaller window ~ poor frequency resolution, smoother PSD since averaged
% over more windows
% Longer window ~ better frequency resolution, noisier PSD since averaging
% over fewer windows

% Welch PSD options
psd_options.n_periods = 24;      % Number of wave periods per window
psd_options.overlap   = 0.5;     % 50% overlap
psd_options.detrend   = 'linear';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE PSD OF EACH SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save options used
spectra.Options = psd_options;

clc; fprintf('Computing PSD...\n\n')
% Loop through farm arrangements
for a = 1:length(farm_arrangements)
    arrangement = farm_arrangements{a};

    % Loop through different farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat("WT60_", farm_spacing, "_AG0");
        fprintf('%s: %s\n', arrangement, caze)
        
        % Get waves available for each case
        waves = fieldnames(full_tracking.(arrangement).(farm_spacing));
    
        % Loop through waves
        for w = 1:length(waves)
            wave = waves{w};
            split_wave = split(wave, '_');
            wavelength_name = split_wave{1};

            % Check if no-wave before getting wave frequency
            if ~strcmp(wavelength_name, 'LM0')
                wave_frequency = forcing_frequencies.(wavelength_name);
            else
                wave_frequency = nan;
            end
    
            % Loop through DOFs
            for d = 1:length(DOFs)
                DOF = DOFs{d};
    
                % Loop through turbines
                for t = 1:length(turbine_catalog.(arrangement).turbines)

                    % Get signal
                    signal = full_tracking.(arrangement).(farm_spacing).(wave)(t).(DOF);
    
                    % Compute power spectral density
                    psd_result = compute_welch_psd(signal, sampling_frequency, wave_frequency, psd_options);

                    % Save results
                    spectra.(arrangement).(farm_spacing).(wave)(t).(DOF) = psd_result;

                end
            end
        end
    end
    fprintf('\n')
end


% CONSIDER ALL NO-WAVE SIGNALS PER DOF TO GET A MEDIAN
% SINGLE SIGNAL

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE PSD OF EACH SIGNAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stack PSDs for all turbines for a given DOF and then take the median at
% each frequency, giving a curve that below which can be considered noise

clear PSD_tmp frequency_tmp
clc; clc; fprintf('Computing Median No-Wave PSD...\n\n')
% Loop through DOFs
for d = 1:length(DOFs)
    DOF = DOFs{d};
    fprintf('%s\n', DOF_names.(DOF))

    % Generate temporay structure to hold all spectra
    PSD_tmp = [];
    frequency_tmp = [];

    % Loop through farm arrangements
    for a = 1:length(farm_arrangements)
        arrangement = farm_arrangements{a};

        % Loop through different farm spacings
        for s = 1:length(farm_spacings)
            farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        
            % Loop through turbines
            for t = 1:length(turbine_catalog.(arrangement).turbines)
    
                % Get PSD
                signal = spectra.(arrangement).(farm_spacing).('LM0_AK00')(t).(DOF).Pxx;
                frequency = spectra.(arrangement).(farm_spacing).('LM0_AK00')(t).(DOF).f;

                % Save to temporary structure
                PSD_tmp = [PSD_tmp; signal.'];
                frequency_tmp = [frequency_tmp; frequency.'];

            end
        end
    end

    % Save to array
    no_wave_median_PSD.(DOF).Pxx = median(PSD_tmp, 1, 'omitnan');
    no_wave_median_PSD.(DOF).f = median(frequency_tmp, 1, 'omitnan');

    % Convert to dB
    no_wave_median_PSD.(DOF).Pxx_dB = pow2db(no_wave_median_PSD.(DOF).Pxx);

end

% Add no-wave median to spectra structure before saving
spectra.NoWaveMedian = no_wave_median_PSD;

% Save STD to matfile for book-keeping and later coupling analysis
fprintf('\n')
save_path = fullfile(tracking_path, 'Spectra', sprintf('OffshoreTracking_PSD_nPeriods_%2.0f_SavitskyGolay_Cutoff_%1.0f.mat', psd_options.n_periods, harmonic_cutoff));
fprintf('Saving Spectral Data...\n')
save(save_path, 'spectra')
fprintf('Done Saving!\n\n')

clear a arrangement c d DOF farm_spacing frequecy frequency_tmp PSD_tmp s signal t
clear s w d t a farm_spacing caze waves wave DOF save_path arrangement signal psd_results
clear split_wave wave_frequency wavelength_name psd_result


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: NO-WAVE SPECTRA
% CENTER TURBINES
% FIXED FARM SPACING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Case to plot
% spacing = 'SX50';
% wave = 'LM0_AK00';
% colors = slanCM(45, length(centers));
% 
% % Plotting
% clc; close all
% figure('color','white')
% t = tiledlayout(2,3);
% sgtitle(sprintf('%s Floating Wind Farm: No-Wave Power Spectra\n$S_x = %1.1fD$', farm_arrangement, spacings_from_string.(spacing)), 'Interpreter', 'latex')
% 
% % Loop through DOF
% for d = 1:length(DOFs)
%     DOF = DOFs{d};
% 
%     % Properly name DOF
%     name = DOF_names.(DOF);
%     symb = DOF_symbs.(DOF);
% 
%     % Plotting
%     h(d) = nexttile;
%     title(sprintf('%s ($%s$)', name, symb), 'interpreter', 'latex', 'fontsize', 14)
%     hold on 
% 
%     % Loop through center turbines
%     for c = 1:length(centers)
% 
%         % Get turbine
%         turbine = centers(c);
% 
%         % Get signals
%         frequency = spectra.(farm_arrangement).(spacing).(wave)(turbine).(DOF).f;
%         power_spectra_dB = spectra.(farm_arrangement).(spacing).(wave)(turbine).(DOF).Pxx_dB;
% 
%         % Plot
%         plot(frequency, power_spectra_dB, ...
%              'color', colors(c,:), 'linewidth', 2, 'HandleVisibility', 'off')
%         scatter(frequency, power_spectra_dB, ...
%              10, 'filled', 'markerfacecolor', colors(c,:), 'HandleVisibility', 'off')
% 
%     end
% 
%     % Legend
%     if d == 1
%         for c = 1:length(centers)
%             plot(nan, nan, 'linewidth', 2, 'Color', colors(c,:), 'DisplayName', sprintf('Row %1.0f', c))
%         end
%         leg = legend('interpreter', 'latex', 'box', 'off');
%         leg.Layout.Tile = 'east';
%     end
% 
%     hold off
% end
% xlabel(t, '$f$', 'Interpreter','latex')
% linkaxes(h, 'x')
% linkaxes(h(1:3), 'y')
% linkaxes(h(4:6), 'y')
% 
% clear d DOF name symb h c turbine t spacing wave row_tag colors power_spectra frequency power_spectra_dB
% clear noise_floor noise_floor_dB tmp

%% Try to find which case + turbine has the bad rotation spectra

% % Case to plot
% spacing = 'SX50';
% wave = 'LM0_AK00';
% colors = slanCM(45, length(turbines));
% 
% % Plotting
% clc; close all
% figure('color','white')
% t = tiledlayout(1,3);
% sgtitle(sprintf('%s Floating Wind Farm: No-Wave Power Spectra\n$S_x = %1.1fD$', farm_arrangement, spacings_from_string.(spacing)), 'Interpreter', 'latex')
% 
% % Loop through DOF
% for d = 1:length(rotations)
%     DOF = rotations{d};
% 
%     % Properly name DOF
%     name = DOF_names.(DOF);
%     symb = DOF_symbs.(DOF);
% 
%     % Plotting
%     h(d) = nexttile;
%     title(sprintf('%s ($%s$)', name, symb), 'interpreter', 'latex', 'fontsize', 14)
%     hold on 
% 
%     % Loop through center turbines
%     for c = 1:length(turbines)
% 
%         % Get turbine
%         turbine = turbines(c);
% 
%         % Get signals
%         frequency = spectra.(farm_arrangement).(spacing).(wave)(turbine).(DOF).f;
%         power_spectra_dB = spectra.(farm_arrangement).(spacing).(wave)(turbine).(DOF).Pxx_dB;
% 
%         % Plot
%         plot(frequency, power_spectra_dB, ...
%              'color', colors(c,:), 'linewidth', 2, 'HandleVisibility', 'off')
%         scatter(frequency, power_spectra_dB, ...
%              10, 'filled', 'markerfacecolor', colors(c,:), 'HandleVisibility', 'off')
% 
%     end
% 
%     % Legend
%     if d == 1
%         for c = 1:length(turbines)
%             plot(nan, nan, 'linewidth', 2, 'Color', colors(c,:), 'DisplayName', sprintf('Turbine %1.0f', c))
%         end
%         leg = legend('interpreter', 'latex', 'box', 'off');
%         leg.Layout.Tile = 'east';
%     end
% 
%     hold off
% end
% xlabel(t, '$f$', 'Interpreter','latex')
% linkaxes(h, 'x')
% linkaxes(h, 'xy')
% 
% clear d DOF name symb h c turbine t spacing wave row_tag colors power_spectra frequency power_spectra_dB
% clear noise_floor noise_floor_dB tmp





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: NO-WAVE SPECTRA + MEDIAN
% CENTER TURBINES
% ALL FARM SPACINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case to plot
wave = 'LM0_AK00';
colors = slanCM(45, length(turbines));

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(2,3);
sgtitle('Inline/Staggered Floating Wind Farm: No-Wave Power Spectra\\nAll Farm Spacings', 'Interpreter', 'latex')

% Loop through DOF
for d = 1:length(DOFs)
    DOF = DOFs{d};
    
    % Properly name DOF
    name = DOF_names.(DOF);
    symb = DOF_symbs.(DOF);
   
    % Plotting
    h(d) = nexttile;
    title(sprintf('%s ($%s$)', name, symb), 'interpreter', 'latex', 'fontsize', 14)
    hold on 

    % Loop through farm arrangements
    for a = 1:length(farm_arrangements)
        arrangement = farm_arrangements{a};

        % Loop through farm spacings
        for s = 1:length(farm_spacings)
            farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
    
            num_turbines = length(turbine_catalog.(arrangement).turbines);
            % Loop through all turbines
            for c = 1:num_turbines
        
                % Get turbine
                turbine = c;
        
                % Get signals
                frequency = spectra.(arrangement).(farm_spacing).(wave)(turbine).(DOF).f;
                power_spectra_dB = spectra.(arrangement).(farm_spacing).(wave)(turbine).(DOF).Pxx_dB;
        
                % Plot
                P = plot(frequency, power_spectra_dB, ...
                         'color', colors(c,:), 'linewidth', 2, 'HandleVisibility', 'off');
                P.Color(4) = 0.05;
                scatter(frequency, power_spectra_dB, ...
                     10, 'filled', 'markerfacecolor', colors(c,:), ...
                     'markerfacealpha', 0.05, 'HandleVisibility', 'off')
                
            end
        end
    end

    % Plot median spectra
    plot(no_wave_median_PSD.(DOF).f, no_wave_median_PSD.(DOF).Pxx_dB, ...
         'color', 'black', 'linewidth', 2)
    scatter(no_wave_median_PSD.(DOF).f, no_wave_median_PSD.(DOF).Pxx_dB, ...
            10, 'filled', 'markerfacecolor', 'black', ...
            'HandleVisibility', 'off')

    
    hold off
    xscale('log')
end
xlabel(t, '$f$ [Hz]', 'Interpreter','latex')
linkaxes(h, 'x')
linkaxes(h(1:3), 'y')
linkaxes(h(4:6), 'y')

clear d DOF name symb h c turbine t spacing wave row_tag colors power_spectra frequency power_spectra_dB
clear noise_floor noise_floor_dB tmp P a arrangement farm_spacing num_turbines s 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: MEDIAN NO-WAVE SPECTRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(2,3);
sgtitle('Inline/Staggered Floating Wind Farm: No-Wave Median Power Spectra', 'Interpreter', 'latex')

% Loop through DOF
for d = 1:length(DOFs)
    DOF = DOFs{d};
    
    % Properly name DOF
    name = DOF_names.(DOF);
    symb = DOF_symbs.(DOF);
   
    % Plotting
    h(d) = nexttile;
    title(sprintf('%s ($%s$)', name, symb), 'interpreter', 'latex', 'fontsize', 14)
    hold on 

    % Get signals
    frequency = no_wave_median_PSD.(DOF).f;
    median_PSD = no_wave_median_PSD.(DOF).Pxx_dB;

    % Last point is ugly, replace with nan
    median_PSD(end) = nan;

    % Plot median spectra
    plot(frequency, median_PSD, ...
         'color', 'black', 'linewidth', 2)
    scatter(frequency, median_PSD, ...
            10, 'filled', 'markerfacecolor', 'black', ...
            'HandleVisibility', 'off')
    hold off

    % Find peaks
    [~, peak_ind] = max(median_PSD, [], 'all', 'omitnan');
    xline(frequency(peak_ind), 'linestyle', '--', ...
          'label', sprintf('$f = %1.2f$ Hz', frequency(peak_ind)), ...
          'LabelVerticalAlignment', 'bottom', ...
          'interpreter', 'latex')
    
    if ismember(DOF, translations)
        ylim([-100, -55])
    elseif ismember(DOF, rotations)
        ylim([-50, 0])
    end
end
xlabel(t, '$f$ [Hz]', 'Interpreter','latex')
linkaxes(h, 'x')
linkaxes(h(1:3), 'y')
linkaxes(h(4:6), 'y')

clear d DOF name symb h c turbine t spacing wave row_tag colors power_spectra frequency power_spectra_dB
clear noise_floor noise_floor_dB tmp peak_ind


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: SPECTRA FOR DIFFERENT WAVES
% CENTER TURBINES
% FREQUENCY NORMALIZED BY WAVE FREQUENCY
% FIXED FARM SPACING
% FIXED DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case to plot
spacing = 'SX50';
DOF = 'yaw_kal';
available_wave_cases = fieldnames(tracking.(spacing));

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(1,length(centers));
sgtitle(sprintf('%s Floating Wind Farm: %s ($%s$) Power Spectra\n$S_x = %1.1fD$', farm_arrangement, name, symb, spacings_from_string.(spacing)), 'Interpreter', 'latex')

% Loop through DOF
for r = 1:length(centers)

    % Get turbine
    turbine = centers(r);
    colors = row_colors.(['Row', num2str(r)]);

    % Plotting
    h(r) = nexttile;
    title(sprintf('Row %1.0f', r), 'interpreter', 'latex', 'fontsize', 14)
    hold on 

    % Loop through steepnesses
    for st = 1:length(wave_steepnesses)
        wave_steepness = wave_steepnesses(st);
        steep = compose('%02d', round(100 * wave_steepness));

        % Loop through wavelengths
        for w = 1:length(wavelengths)
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            % Check that wave case is available in tracking data
            if ismember(wave, available_wave_cases)

                % Apply legend to only one wave steepness
                if st == 3
                    vis = 'on';
                else
                    vis = 'off';
                end
           
                % Get data
                normalized_frequency = spectra.(farm_arrangement).(spacing).(wave)(turbine).(DOF).f_norm;
                power_spectra_dB = spectra.(farm_arrangement).(spacing).(wave)(turbine).(DOF).Pxx_dB;
        
                % Plot lines
                P = plot(normalized_frequency, power_spectra_dB, ...
                         'color', colors(w,:), 'linewidth', 2, ...
                         'HandleVisibility', vis, ...
                         'Displayname', sprintf('$\\lambda = %1.0f$', wavelengths(w)));
                P.Color(4) = steepness_alpha(st);

                % Plot points
                scatter(normalized_frequency, power_spectra_dB, ...
                     10, 'filled', 'markerfacecolor', colors(w,:), 'markerfacealpha', steepness_alpha(st), ...
                     'HandleVisibility', 'off')
            end
        end
    end
    hold off
    legend('interpreter', 'latex', 'box', 'off', 'location', 'northeast', 'fontsize', 8)
    xlim([0, 3])

    if ismember(DOF, rotations)
        ylim([-60, 40])
    else
        ylim([-130, -20])
    end

    % Add markers at harmonics
    for i = 0.5:0.5:3
        xline(i, 'linestyle', '--', 'HandleVisibility', 'off', 'alpha', 0.5)
    end

end
xlabel(t, '$f / f_{wave}$', 'Interpreter','latex')
ylabel(t, sprintf('%s ($%s$) PSD dB', name, symb), 'Interpreter','latex')
linkaxes(h, 'x')
linkaxes(h, 'xy')


clear d DOF name symb h c turbine t spacing wave row_tag colors available_wave_cases flow_advection_frequecy i P r
clear st steep vis w wave_steepness


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: SPECTRA FOR DIFFERENT WAVES, CENTER TURBINES
% FREQUENCY NOT-NORMALIZED
% COMPARING TO NO-WAVE MEDIAN
% FIXED FARM SPACING
% FIXED DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case to plot
spacing = 'SX50';
DOF = 'yaw_kal';
available_wave_cases = fieldnames(tracking.(spacing));

% Set y-limits based on DOF
if ismember(DOF, translations)
    ylims = [-120, -20];
elseif ismember(DOF, rotations)
    ylims = [-60, 40];
end

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(1,length(centers));
sgtitle(sprintf('%s Floating Wind Farm: %s ($%s$) Power Spectra\n$S_x = %1.1fD$', farm_arrangement, name, symb, spacings_from_string.(spacing)), 'Interpreter', 'latex')

% Loop through DOF
for r = 1:length(centers)

    % Get turbine
    turbine = centers(r);
    colors = row_colors.(['Row', num2str(r)]);

    % Plotting
    h(r) = nexttile;
    title(sprintf('Row %1.0f', r), 'interpreter', 'latex', 'fontsize', 14)
    hold on 

    % Loop through steepnesses
    for st = 1:length(wave_steepnesses)
    % for st = 3
        wave_steepness = wave_steepnesses(st);
        steep = compose('%02d', round(100 * wave_steepness));

        % Loop through wavelengths
        for w = 1:length(wavelengths)
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            % Check that wave case is available in tracking data
            if ismember(wave, available_wave_cases)

                % Apply legend to only one wave steepness
                if st == 3
                    vis = 'on';
                else
                    vis = 'off';
                end
           
                % Get data
                normalized_frequency = spectra.(farm_arrangement).(spacing).(wave)(turbine).(DOF).f;
                power_spectra_dB = spectra.(farm_arrangement).(spacing).(wave)(turbine).(DOF).Pxx_dB;

                % No wave data
                no_wave_frequency = no_wave_median_PSD.(DOF).f;
                median_PSD = no_wave_median_PSD.(DOF).Pxx_dB;

                if st == 3 && w == 1
                    no_wave_vis = 'on';
                else
                    no_wave_vis = 'off';
                end

                % Plot no-wave median
                plot(no_wave_frequency, median_PSD, ...
                    'color', 'black', 'linewidth', 2, ...
                    'displayname', 'No-Wave', ...
                    'HandleVisibility', no_wave_vis);

                scatter(no_wave_frequency, median_PSD, ...
                     10, 'filled', 'markerfacecolor', 'black', ...
                     'HandleVisibility', 'off')



                % Shade underneath no-wave
                if st == 3 && w == 1
                    % Construct patch coordinates
                    ylim(ylims)
                    ybase = ylims(1);
                    x_patch = [no_wave_frequency(:); flipud(no_wave_frequency(:))];
                    y_patch = [median_PSD(:); ybase*ones(size(median_PSD(:)))];
                    
                    % Draw patch
                    p = patch(x_patch, y_patch, 'k', ...
                              'FaceAlpha', 0.5, ...
                              'EdgeColor', 'none', ...
                              'HandleVisibility', 'off');
                end


        
                % Plot lines
                P = plot(normalized_frequency, power_spectra_dB, ...
                         'color', colors(w,:), 'linewidth', 2, ...
                         'HandleVisibility', vis, ...
                         'Displayname', sprintf('$\\lambda = %1.0f$', wavelengths(w)));
                P.Color(4) = steepness_alpha(st);

                % Plot points
                scatter(normalized_frequency, power_spectra_dB, ...
                     10, 'filled', 'markerfacecolor', colors(w,:), 'markerfacealpha', steepness_alpha(st), ...
                     'HandleVisibility', 'off')

                if st == 3
                    % Add a marker at wave forcing frequency
                     xline(forcing_frequencies.(['LM', num2str(wavelengths(w))]), ...
                           'linestyle', '--', 'color', colors(w,:), ...
                           'linewidth', 1, 'HandleVisibility', 'off')
                end
            end
        end
    end
    hold off
    legend('interpreter', 'latex', 'box', 'off', 'location', 'northeast', 'fontsize', 8)

    % Set frequecy range
    xlim([0, 6])
end
xlabel(t, '$f$ [Hz]', 'Interpreter','latex')
ylabel(t, sprintf('%s ($%s$) PSD dB', name, symb), 'Interpreter','latex')
linkaxes(h, 'xy')


clear d DOF name symb h c turbine t spacing wave row_tag colors available_wave_cases flow_advection_frequecy i P r
clear st steep vis w wave_steepness








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psd_out] = compute_welch_psd(signal, fs, f_wave, options)
% COMPUTE_WELCH_PSD Compute power spectral density using Welch's method
%
%   [psd_out] = compute_welch_psd(signal, fs, f_wave)
%   [psd_out] = compute_welch_psd(signal, fs, f_wave, options)
%
%   Computes the PSD of a signal using Welch's method with parameters
%   chosen to give good frequency resolution around the wave frequency.
%
%   INPUTS:
%       signal  - Time series data (column vector)
%       fs      - Sampling frequency (Hz)
%       f_wave  - Wave forcing frequency (Hz), used for setting resolution
%                 Use NaN or [] for no-wave cases (will use default window)
%       options - (optional) struct with fields:
%           .n_periods       - Number of wave periods per window (default: 8)
%           .overlap         - Fractional overlap between windows (default: 0.5)
%           .detrend         - Detrend method: 'linear', 'constant', or 'none' (default: 'linear')
%           .default_window  - Window duration in seconds for no-wave cases (default: 4)
%
%   OUTPUTS:
%       psd_out - Struct containing:
%           .f          - Frequency vector (Hz)
%           .f_norm     - Normalized frequency vector (f / f_wave), empty for no-wave cases
%           .Pxx        - Power spectral density (units^2/Hz)
%           .Pxx_dB     - PSD in decibels: 10*log10(Pxx)
%           .window_len - Window length used (samples)
%           .n_windows  - Approximate number of windows averaged
%           .df         - Frequency resolution (Hz)
%           .fs         - Sampling frequency (stored for reference)
%           .f_wave     - Wave frequency (stored for reference, NaN if no wave)
%
%   Zein Sadek, January 2025

    % Default options
    if nargin < 4
        options = struct();
    end
    if ~isfield(options, 'n_periods')
        options.n_periods = 8;  % 8 wave periods per window
    end
    if ~isfield(options, 'overlap')
        options.overlap = 0.5;  % 50% overlap
    end
    if ~isfield(options, 'detrend')
        options.detrend = 'linear';
    end
    if ~isfield(options, 'default_window')
        options.default_window = 4;  % 4 seconds for no-wave cases
    end

    % Ensure column vector
    signal = signal(:);
    
    % Remove NaNs if present (take longest continuous segment)
    if any(isnan(signal))
        warning('Signal contains NaNs - using longest continuous segment');
        nan_idx = isnan(signal);
        segments = diff([0; nan_idx; 0]);
        seg_starts = find(segments == -1);
        seg_ends = find(segments == 1) - 1;
        seg_lengths = seg_ends - seg_starts + 1;
        [~, longest] = max(seg_lengths);
        signal = signal(seg_starts(longest):seg_ends(longest));
    end

    % Detrend
    switch options.detrend
        case 'linear'
            signal = detrend(signal, 'linear');
        case 'constant'
            signal = detrend(signal, 'constant');
        case 'none'
            % Do nothing
    end

    % Determine window length
    % If no wave frequency provided, use default window duration
    if isempty(f_wave) || isnan(f_wave)
        window_time = options.default_window;
        f_wave = NaN;  % Store as NaN for reference
    else
        % Window length: capture n_periods of the wave frequency
        T_wave = 1 / f_wave;                       % Wave period (s)
        window_time = options.n_periods * T_wave;  % Window duration (s)
    end
    window_len = round(window_time * fs);          % Window length (samples)
    
    % Ensure window length doesn't exceed signal length
    N = length(signal);
    if window_len > N
        warning('Window length exceeds signal length. Using N/4 instead.');
        window_len = floor(N / 4);
    end
    
    % Overlap in samples
    n_overlap = round(options.overlap * window_len);
    
    % Compute PSD using Welch's method with Hamming window
    [Pxx, f] = pwelch(signal, hamming(window_len), n_overlap, [], fs);
    
    % Approximate number of windows (for reference)
    step = window_len - n_overlap;
    n_windows = floor((N - n_overlap) / step);
    
    % Frequency resolution
    df = f(2) - f(1);
    
    % Package outputs
    psd_out.f          = f;
    psd_out.Pxx        = Pxx;
    psd_out.Pxx_dB     = 10 * log10(Pxx + eps);  % eps to avoid log(0)
    psd_out.window_len = window_len;
    psd_out.n_windows  = n_windows;
    psd_out.df         = df;
    psd_out.fs         = fs;
    psd_out.f_wave     = f_wave;
    
    % Normalized frequency (f / f_wave) - only for wave cases
    if ~isnan(f_wave)
        psd_out.f_norm = f / f_wave;
    else
        psd_out.f_norm = [];  % Empty for no-wave cases
    end

end