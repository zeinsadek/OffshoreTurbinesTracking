%% Looking at band-limited RMS/STD of motion (CORRECTED VERSION)
% Zein Sadek

% Key changes from original:
% 1. Uses FFT-based band-limited RMS (no filter roll-off losses)
% 2. Bands are defined to be perfectly contiguous: 
% [0, 0.5*fw], [0.5*fw, 1.5*fw], [1.5*fw, Fs/2]
% 3. Includes verification that sum of squared band RMS equals total variance

% Saved matfile is technically decomposed standard deviation (having same
% units as signal) and need to be converted into variance before being
% combined: full^2 = LF^2 + WF^2 + HF^2, hence the same
% bandPartitionedDEVIATIONS

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

% Load full STDs
deviations_file = sprintf('OffshoreTracking_FullSTD_SavitskyGolay_Cutoff_%1.0f.mat', harmonic_cutoff);
full_deviations = load(fullfile(tracking_path, "StandardDeviations", deviations_file));
full_deviations = full_deviations.deviations;

% Load data for specific layout
deviations = full_deviations.(farm_arrangement);

clear offshore_path tracking_file deviations_file


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

% Frequency component names
frequency_names.LF = "Low Frequency"; 
frequency_names.WF = "Wave Frequency"; 
frequency_names.HF = "High Frequency"; 


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
% COMPUTE BAND-LIMITED RMS OF EACH DOF (FFT METHOD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KEY CHANGES:
% 1. Using FFT-based method for sharp frequency cutoffs
% 2. Bands are perfectly contiguous: [0, 0.5*fw], [0.5*fw, 1.5*fw], [1.5*fw, Fs/2]
% 3. This ensures sum of squared RMS values = total variance

sampling_frequency = 30;

clc; fprintf('Computing band-limited RMS using FFT method...\n\n');
% Loop through farm arrangements
for a = 1:length(farm_arrangements)
    arrangement = farm_arrangements{a};

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
    
            % Skip no-wave case
            if ~strcmp(wavelength_name, 'LM0')
                fw = forcing_frequencies.(wavelength_name);
    
                % Define contiguous bands spanning [0, Fs/2]
                % No gaps, no overlaps
                bands.LF = [0,       0.5*fw];                 % DC to half wave frequency
                bands.WV = [0.5*fw,  1.5*fw];                 % Wave frequency band
                bands.HF = [1.5*fw,  sampling_frequency/2];   % High frequency to Nyquist
        
                % Loop through DOFs
                for d = 1:length(DOFs)
                    DOF = DOFs{d};
        
                    % Loop through turbines
                    for t = 1:length(turbine_catalog.(arrangement).turbines)
                        % Signal
                        q = full_tracking.(arrangement).(farm_spacing).(wave)(t).(DOF);
        
                        % Compute low frequency portion (includes DC)
                        rLF = bandlimited_rms_fft(q, sampling_frequency, bands.LF(1), bands.LF(2), 'constant');
                        
                        % Compute wave frequency portion
                        rWV = bandlimited_rms_fft(q, sampling_frequency, bands.WV(1), bands.WV(2), 'constant');
        
                        % Compute high frequency portion (up to Nyquist)
                        rHF = bandlimited_rms_fft(q, sampling_frequency, bands.HF(1), bands.HF(2), 'constant');
    
                        % Save
                        bandPartitionedDeviations.(arrangement).(farm_spacing).(wave)(t).(DOF).LF = rLF;
                        bandPartitionedDeviations.(arrangement).(farm_spacing).(wave)(t).(DOF).WF = rWV;
                        bandPartitionedDeviations.(arrangement).(farm_spacing).(wave)(t).(DOF).HF = rHF;
                        bandPartitionedDeviations.(arrangement).(farm_spacing).(wave)(t).(DOF).Full = full_deviations.(arrangement).(farm_spacing).(wave)(t).(DOF);
       
                    end
                end
            end
        end
    end
    fprintf('\n')
end
fprintf('Done!\n');

clear a arrangement s farm_spacing caze waves w wave split_wave wavelength_name fw
clear bands d DOF t q rLF rWV rHF


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK ALL CASES - Find worst closure errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ratios = [];
worst_cases = struct('ratio', [], 'spacing', {}, 'wave', {}, 'turbine', [], 'DOF', {});

clc; fprintf('=== Checking closure errors for all cases ===\n\n');
% Loop through farm arrangements
for a = 1:length(farm_arrangements)
    arrangement = farm_arrangements{a};

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        waves = fieldnames(full_tracking.(arrangement).(farm_spacing));
    
        % Loop through waves
        for w = 1:length(waves)
            wave = waves{w};
            split_wave = split(wave, '_');
            wavelength_name = split_wave{1};
    
            % Skip no-wave cases
            if ~strcmp(wavelength_name, 'LM0') && isfield(bandPartitionedDeviations.(arrangement).(farm_spacing), wave)

                % Loop through DOFs
                for d = 1:length(DOFs)
                    DOF = DOFs{d};

                    % Loop through turbines
                    for t = 1:length(turbine_catalog.(arrangement).turbines)
                        % Get components
                        full = full_deviations.(arrangement).(farm_spacing).(wave)(t).(DOF);
                        LF = bandPartitionedDeviations.(arrangement).(farm_spacing).(wave)(t).(DOF).LF;
                        WF = bandPartitionedDeviations.(arrangement).(farm_spacing).(wave)(t).(DOF).WF;
                        HF = bandPartitionedDeviations.(arrangement).(farm_spacing).(wave)(t).(DOF).HF;
                        
                        % Check summation
                        ss = LF^2 + WF^2 + HF^2;
                        ratio = ss / full^2;
                        ratios = [ratios; ratio];
                        
                        % Track worst cases
                        if abs(ratio - 1) > 0.01  % More than 1% error
                            idx = length(worst_cases) + 1;
                            worst_cases(idx).ratio = ratio;
                            worst_cases(idx).spacing = farm_spacing;
                            worst_cases(idx).wave = wave;
                            worst_cases(idx).turbine = t;
                            worst_cases(idx).DOF = DOF;
                        end
                    end
                end
            end
        end
    end
end

% Display statistics
fprintf('Closure ratio statistics:\n');
fprintf('  Min:    %.6f\n', min(ratios));
fprintf('  Max:    %.6f\n', max(ratios));
fprintf('  Mean:   %.6f\n', mean(ratios));
fprintf('  Median: %.6f\n', median(ratios));
fprintf('  Std:    %.6f\n', std(ratios));
fprintf('\nCases with >1%% error: %d out of %d\n\n', length(worst_cases), length(ratios));

if ~isempty(worst_cases)
    fprintf('\nWorst cases:\n');
    [~, sortIdx] = sort(abs([worst_cases.ratio] - 1), 'descend');
    for i = 1:min(5, length(sortIdx))
        wc = worst_cases(sortIdx(i));
        fprintf('  %s, %s, T%d, %s: ratio = %.4f\n', ...
            wc.spacing, wc.wave, wc.turbine, wc.DOF, wc.ratio);
    end
end

% Save matfile
fprintf('\n')
save_path = fullfile(tracking_path, 'StandardDeviations', sprintf('OffshoreTracking_BandPartitionedSTD_SavitskyGolay_Cutoff_%1.0f.mat', harmonic_cutoff));
fprintf('Saving Band-Partitioned Standard Deviations...\n')
save(save_path, 'bandPartitionedDeviations')
fprintf('Done Saving!\n\n')

% Only consider deviations from specific arrangement
bandPartitionedDeviations = bandPartitionedDeviations.(farm_arrangement);

clear a arrangement s farm_spacing waves w wave split_wave wavelength_name d DOF t full
clear LF WF HF ss ratio ratios worst_cases wc sortIdx


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: NORMALIZED BAND-PARTITIONED VARIANCE 
% AGAINST WAVELENGTH
% LOOPED OVER ALL DOF
% CONSIDERING A SINGLE TURBINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the variance associated with the low, wave, or high 
% frequency bands normalized by the total variance
% Represents the percentage of the energy contained in this band, and varys
% from 0 to 1

% Which turbine to plot
turbine = centers(3);
frequency_component = 'HF';
row = turbineRow(turbine, farm_arrangement);
frequency_name = frequency_names.(frequency_component);
colors = row_colors.(sprintf('Row%1.0f', row));

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(2,3);
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f %s Wave Score: $\\sigma_{%s}^2 / \\sigma^2$', farm_arrangement, row, frequency_name, frequency_component), 'Interpreter', 'latex')

% Loop over DOFs
for d = 1:length(DOFs)
    DOF = DOFs{d};

    % Properly name DOF
    name = DOF_names.(DOF);
    symb = DOF_symbs.(DOF);    

    % Plotting
    h(d) = nexttile;
    title(sprintf('%s: $(%s)$', name, symb), 'interpreter', 'latex', 'fontsize', 14)
    hold on 

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(bandPartitionedDeviations.(farm_spacing));

        % Loop through wave steepnesses
        for st = 1:length(wave_steepnesses)
            steep = compose('%02d', round(100 * wave_steepnesses(st)));
     
            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)
    
                    total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
                    wave_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).(frequency_component)^2;
                    wave_score = wave_band_variance ./ total_variance;
    
                    scatter(wavelengths(w), wave_score, marker_size, 'filled', ...
                            'MarkerFaceColor', colors(s,:), ...
                            'MarkerFaceAlpha', steepness_alpha(st), ...
                            'Marker', wave_shapes{w}, ...
                            'HandleVisibility', 'off')

                    tmpX(w) = wavelengths(w);
                    tmpY(w) = wave_score;
                end
            end
            % Connect the dots
            P = plot(tmpX, tmpY, 'color', colors(s,:), 'linewidth', 1, 'HandleVisibility', 'off');
            P.Color(4) = steepness_alpha(st);
        end
    end

    % Legend
    if d == 1
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
xlabel(t, '$\lambda / D$', 'Interpreter','latex')
ylabel(t, sprintf('$\\sigma_{%s}^2 / \\sigma^2$', frequency_component), 'interpreter', 'latex')
linkaxes(h, 'xy')
xlim([1.5, 5.5])
ylim([0, 1])

clear frequency_component row frequency_name colors t d DOF name symb h s farm_spacing
clear available_wave_cases st wave_steepness steep w wave harmonic_ratio total_variance
clear wave_band_variance wave_score turbine leg tmpX tmpY


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: NORMALIZED BAND-PARTITIONED VARIANCE 
% AGAINST HARMONIC RATIO
% LOOPED OVER ALL DOF
% CONSIDERING A SINGLE TURBINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the variance associated with the low, wave, or high 
% frequency bands normalized by the total variance
% Represents the percentage of the energy contained in this band, and varys
% from 0 to 1

% % Which turbine to plot
% turbine = 5;
% frequency_component = 'LF';
% row = turbineRow(turbine, farm_arrangement);
% frequency_name = frequency_names.(frequency_component);
% colors = row_colors.(sprintf('Row%1.0f', row));
% 
% % Plotting
% clc; close all
% figure('color','white')
% t = tiledlayout(2,3);
% sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f %s Wave Score: $\\sigma_{%s}^2 / \\sigma^2$', farm_arrangement, row, frequency_name, frequency_component), 'Interpreter', 'latex')
% 
% % Loop over DOFs
% for d = 1:length(DOFs)
%     DOF = DOFs{d};
% 
%     % Properly name DOF
%     name = DOF_names.(DOF);
%     symb = DOF_symbs.(DOF);    
% 
%     % Plotting
%     h(d) = nexttile;
%     title(sprintf('%s: $(%s)$', name, symb), 'interpreter', 'latex', 'fontsize', 14)
%     hold on 
% 
%     % Loop through farm spacings
%     for s = 1:length(farm_spacings)
%         farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%         available_wave_cases = fieldnames(bandPartitionedDeviations.(farm_spacing));
% 
%         % Loop through wave steepnesses
%         for st = 1:length(wave_steepnesses)
%             steep = compose('%02d', round(100 * wave_steepnesses(st)));
% 
%             % Loop through wavelengths
%             for w = 1:length(wavelengths)
%                 wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
% 
%                 % Check that wave case is available in tracking data
%                 if ismember(wave, available_wave_cases)
%                     harmonic_ratio = farm_spacings(s) / wavelengths(w);
% 
%                     total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
%                     wave_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).(frequency_component)^2;
%                     wave_score = wave_band_variance ./ total_variance;
% 
%                     scatter(harmonic_ratio, wave_score, marker_size, 'filled', ...
%                             'MarkerFaceColor', colors(s,:), ...
%                             'MarkerFaceAlpha', steepness_alpha(st), ...
%                             'Marker', wave_shapes{w}, ...
%                             'HandleVisibility', 'off')
%                 end
%             end
%         end
%     end
% 
%     % Legend
%     if d == 1
%         % Legend for color
%         for s = 1:length(farm_spacings)
%             plot(nan, nan, 'Color', colors(s,:), 'linewidth', 3, ...
%                 'Displayname', sprintf('$S_x = %1.1fD', farm_spacings(s)), 'HandleVisibility', 'on')
%         end
% 
%         % White space
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
% 
%         % Legend for marker shape
%         for w = 1:length(wavelengths)
%             scatter(nan, nan, marker_size, wave_shapes{w}, 'black', 'filled', 'HandleVisibility', 'on', ...
%                     'DisplayName', sprintf('$\\lambda = %1.0fD$', wavelengths(w)))
%         end
% 
%         % White space
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
% 
%         % Legend for marker alpha
%         for st = 1:length(wave_steepnesses)
%             scatter(nan, nan, marker_size, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
%                     'markerfacealpha', steepness_alpha(st), ...
%                     'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
%         end
% 
%         leg = legend('interpreter', 'latex', 'box', 'off');
%         leg.Layout.Tile = 'east';
%     end
%     hold off
% end
% xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
% ylabel(t, sprintf('$\\sigma_{%s}^2 / \\sigma^2$', frequency_component), 'interpreter', 'latex')
% linkaxes(h, 'xy')
% xlim([0.5, 2.6])
% ylim([0, 1])
% 
% clear frequency_component row frequency_name colors t d DOF name symb h s farm_spacing
% clear available_wave_cases st wave_steepness steep w wave harmonic_ratio total_variance
% clear wave_band_variance wave_score turbine leg


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: BALANCE BETWWEEN LOW AND WAVE COMPONENTS 
% BAND-PARTITIONED VARIANCE AGAINST WAVELENGTH
% LOOPED OVER ALL DOF
% CONSIDERING A SINGLE TURBINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the variance associated with the low, wave, or high 
% frequency bands normalized by the total variance
% Represents the balance of the energy contained in the low and wave bands
% Varys from -1 to 1: -1 ~ low frequency, +1 ~ wave frequency

% Which turbine to plot
turbine = centers(3);
row = turbineRow(turbine, farm_arrangement);
colors = row_colors.(sprintf('Row%1.0f', row));

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(2,3);
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f Low/Wave Frequency Variance Balance', farm_arrangement, row), 'Interpreter', 'latex')

% Loop through DOFs
for d = 1:length(DOFs)
    DOF = DOFs{d};

    % Properly name DOF
    name = DOF_names.(DOF);
    symb = DOF_symbs.(DOF);    
       
    % Plotting
    clc; clear tmp
    h(d) = nexttile;
    title(sprintf('%s: $(%s)$', name, symb), 'interpreter', 'latex', 'fontsize', 14)
    hold on 

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(bandPartitionedDeviations.(farm_spacing));

        % Loop through wave steepnesses
        for st = 1:length(wave_steepnesses)
            steep = compose('%02d', round(100 * wave_steepnesses(st)));
 
            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)
    
                    % Get data
                    total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
                    low_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).LF^2;
                    wave_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).WF^2;
    
                    % Normalize
                    low_norm = low_band_variance / total_variance;
                    wave_norm = wave_band_variance / total_variance;
    
                    % Compute score
                    balance = (wave_norm - low_norm) / (wave_norm + low_norm);
    
                    scatter(wavelengths(w), balance, marker_size, 'filled', ...
                            'MarkerFaceColor', colors(s,:), ...
                            'MarkerFaceAlpha', steepness_alpha(st), ...
                            'Marker', wave_shapes{w}, ...
                            'HandleVisibility', 'off')

                    tmpX(w) = wavelengths(w);
                    tmpY(w) = balance;
                end
            end
            % Connect the dots
            P = plot(tmpX, tmpY, 'color', colors(s,:), 'linewidth', 1, 'HandleVisibility', 'off');
            P.Color(4) = steepness_alpha(st);
        end
    end

    % Legend
    if d == 1
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
xlabel(t, '$\lambda / D$', 'Interpreter','latex')
ylabel(t, '$(P_{WF} - P_{LF}) / (P_{WF} + P_{LF})$', 'interpreter', 'latex')
linkaxes(h, 'xy')
xlim([1.5, 5.5])
ylim([-1, 1])


clear frequency_component row frequency_name colors t d DOF name symb h s farm_spacing
clear available_wave_cases st wave_steepness steep w wave harmonic_ratio total_variance 
clear wave_band_variance wave_score turbine total_variance low_band_variance wave_band_variance
clear low_norm wave_norm balance leg


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: BALANCE BETWWEEN LOW AND WAVE COMPONENTS 
% BAND-PARTITIONED VARIANCE AGAINST HARMONIC RATIO
% LOOPED OVER ALL DOF
% CONSIDERING A SINGLE TURBINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the variance associated with the low, wave, or high 
% frequency bands normalized by the total variance
% Represents the balance of the energy contained in the low and wave bands
% Varys from -1 to 1: -1 ~ low frequency, +1 ~ wave frequency

% % Which turbine to plot
% turbine = 5;
% row = turbineRow(turbine, farm_arrangement);
% colors = row_colors.(sprintf('Row%1.0f', row));
% 
% % Plotting
% clc; close all
% figure('color','white')
% t = tiledlayout(2,3);
% sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f Low/Wave Frequency Variance Balance', farm_arrangement, row), 'Interpreter', 'latex')
% 
% % Loop through DOFs
% for d = 1:length(DOFs)
%     DOF = DOFs{d};
% 
%     % Properly name DOF
%     name = DOF_names.(DOF);
%     symb = DOF_symbs.(DOF);    
% 
%     % Plotting
%     clc; clear tmp
%     h(d) = nexttile;
%     title(sprintf('%s: $(%s)$', name, symb), 'interpreter', 'latex', 'fontsize', 14)
%     hold on 
% 
%     % Loop through farm spacings
%     for s = 1:length(farm_spacings)
%         farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%         available_wave_cases = fieldnames(bandPartitionedDeviations.(farm_spacing));
% 
%         % Loop through wave steepnesses
%         for st = 1:length(wave_steepnesses)
%             steep = compose('%02d', round(100 * wave_steepnesses(st)));
% 
%             % Loop through wavelengths
%             for w = 1:length(wavelengths)
%                 wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
% 
%                 % Check that wave case is available in tracking data
%                 if ismember(wave, available_wave_cases)
%                     harmonic_ratio = farm_spacings(s) / wavelengths(w);
% 
%                     % Get data
%                     total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
%                     low_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).LF^2;
%                     wave_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).WF^2;
% 
%                     % Normalize
%                     low_norm = low_band_variance / total_variance;
%                     wave_norm = wave_band_variance / total_variance;
% 
%                     % Compute score
%                     balance = (wave_norm - low_norm) / (wave_norm + low_norm);
% 
%                     scatter(harmonic_ratio, balance, marker_size, 'filled', ...
%                             'MarkerFaceColor', colors(s,:), ...
%                             'MarkerFaceAlpha', steepness_alpha(st), ...
%                             'Marker', wave_shapes{w}, ...
%                             'HandleVisibility', 'off')
%                 end
%             end
%         end
%     end
% 
%     % Legend
%     if d == 1
%         % Legend for color
%         for s = 1:length(farm_spacings)
%             plot(nan, nan, 'Color', colors(s,:), 'linewidth', 3, ...
%                 'Displayname', sprintf('$S_x = %1.1fD', farm_spacings(s)), 'HandleVisibility', 'on')
%         end
% 
%         % White space
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
% 
%         % Legend for marker shape
%         for w = 1:length(wavelengths)
%             scatter(nan, nan, marker_size, wave_shapes{w}, 'black', 'filled', 'HandleVisibility', 'on', ...
%                     'DisplayName', sprintf('$\\lambda = %1.0fD$', wavelengths(w)))
%         end
% 
%         % White space
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
% 
%         % Legend for marker alpha
%         for st = 1:length(wave_steepnesses)
%             scatter(nan, nan, marker_size, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
%                     'markerfacealpha', steepness_alpha(st), ...
%                     'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
%         end
% 
%         leg = legend('interpreter', 'latex', 'box', 'off');
%         leg.Layout.Tile = 'east';
%     end
%     hold off
% 
%     % Set aspect ratio equal since both values are unitless
%     % axis equal
% end
% xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
% ylabel(t, '$(P_{WF} - P_{LF}) / (P_{WF} + P_{LF})$', 'interpreter', 'latex')
% linkaxes(h, 'xy')
% xlim([0.5, 2.6])
% ylim([-1, 1])
% 
% 
% clear frequency_component row frequency_name colors t d DOF name symb h s farm_spacing
% clear available_wave_cases st wave_steepness steep w wave harmonic_ratio total_variance 
% clear wave_band_variance wave_score turbine total_variance low_band_variance wave_band_variance
% clear low_norm wave_norm balance leg


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: BALANCE BETWWEEN LOW AND WAVE COMPONENTS 
% BAND-PARTITIONED VARIANCE AGAINST WAVELENGTH
% LOOPED OVER ALL CENTER TURBINES
% CONSIDERING A DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the variance associated with the low, wave, or high 
% frequency bands normalized by the total variance
% Represents the balance of the energy contained in the low and wave bands
% Varys from -1 to 1: -1 ~ low frequency, +1 ~ wave frequency

% Which DOF to consider
DOF = 'pitch_kal';

% Properly name DOF
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(1,length(centers) + 1);
sgtitle(sprintf('%s Floating Wind Farm: %s ($%s$) Low/Wave Frequency Variance Balance', farm_arrangement, name, symb), 'Interpreter', 'latex')

% Loop through center turbines
for c = 1:length(centers)

    % Get turbine and corresponding color
    turbine = centers(c);
    row = turbineRow(turbine, farm_arrangement);
    colors = row_colors.(sprintf('Row%1.0f', row));

    % Plotting
    h(c) = nexttile;
    title(sprintf('Row %1.0f', c))
    hold on 

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(bandPartitionedDeviations.(farm_spacing));

        % Loop through wave steepnesses
        for st = 1:length(wave_steepnesses)
            steep = compose('%02d', round(100 * wave_steepnesses(st)));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)
    
                    % Get data
                    total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
                    low_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).LF^2;
                    wave_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).WF^2;
    
                    % Normalize
                    low_norm = low_band_variance / total_variance;
                    wave_norm = wave_band_variance / total_variance;
    
                    % Compute score
                    wave_score = (wave_norm - low_norm) / (wave_norm + low_norm);
    
                    scatter(wavelengths(w), wave_score, marker_size, 'filled', ...
                            'MarkerFaceColor', colors(s,:), ...
                            'MarkerFaceAlpha', steepness_alpha(st), ...
                            'Marker', wave_shapes{w}, ...
                            'HandleVisibility', 'off')

                    tmpX(w) = wavelengths(w);
                    tmpY(w) = wave_score;
                end
            end
            % Connect the dots
            P = plot(tmpX, tmpY, 'color', colors(s,:), 'linewidth', 1, 'HandleVisibility', 'off');
            P.Color(4) = steepness_alpha(st);
        end
    end

    % Legend for each plot showing colors
    for ww = 1:length(wavelengths)
        plot(nan, nan, 'Color', colors(ww,:), 'linewidth', 3, ...
            'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(ww)), 'HandleVisibility', 'on')
    end
    legend('Interpreter', 'latex', 'box', 'off', 'location', 'best', 'fontsize', 8)

    % Horizontal line at zero
    yline(0, 'linestyle', '--', 'linewidth', 1, 'HandleVisibility', 'off', 'Alpha', 0.25)
    yticks(-1:0.5:1)
    hold off
    xticks(2:1:5)

    % Set aspect ratio equal since both values are unitless
    % axis equal
end

% Generate outside legend for marker shape and transparency
axLeg = nexttile;
hold on; axis off

% Proxies for marker shapes (spacing)
for ww = 1:numel(wavelengths)
    scatter(axLeg, nan, nan, marker_size, wave_shapes{ww}, 'k', 'filled', ...
        'DisplayName', sprintf('$\\lambda = %1.0fD$', wavelengths(ww)));
end

% Separator
plot(axLeg, nan, nan, 'w', 'DisplayName', ' ');

% Proxies for alpha (steepness)
for stt = 1:numel(wave_steepnesses)
    scatter(axLeg, nan, nan, marker_size, 'o', 'k', 'filled', ...
        'MarkerFaceAlpha', steepness_alpha(stt), ...
        'DisplayName', sprintf('$ak = %1.2f$', wave_steepnesses(stt)));
end

legend('Interpreter', 'latex', 'Box', 'off', 'FontSize', 10, 'Location', 'west');
hold off
  

% Plot axes labels
xlabel(t, '$\lambda / D$', 'Interpreter','latex')
ylabel(t, '$(P_{WF} - P_{LF}) / (P_{WF} + P_{LF})$', 'interpreter', 'latex')
linkaxes(h, 'xy')
xlim(h(1:length(centers)), [1.5, 5.5])
ylim(h(1:length(centers)), [-1, 1])

clear frequency_component row frequency_name colors t d DOF name symb h s farm_spacing
clear available_wave_cases st wave_steepness steep w wave harmonic_ratio total_variance 
clear wave_band_variance wave_score turbine total_variance low_band_variance wave_band_variance
clear low_norm wave_norm balance leg alpha_proxy axLeg c H legGlobal ss stt ww


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: BALANCE BETWWEEN LOW AND WAVE COMPONENTS 
% BAND-PARTITIONED VARIANCE AGAINST HARMONIC RATIO
% LOOPED OVER ALL CENTER TURBINES
% CONSIDERING A DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the variance associated with the low, wave, or high 
% frequency bands normalized by the total variance
% Represents the balance of the energy contained in the low and wave bands
% Varys from -1 to 1: -1 ~ low frequency, +1 ~ wave frequency

% % Which DOF to consider
% DOF = 'pitch_kal';
% 
% % Properly name DOF
% name = DOF_names.(DOF);
% symb = DOF_symbs.(DOF);
% 
% % Plotting
% clc; close all
% figure('color','white')
% t = tiledlayout(1,length(centers) + 1);
% sgtitle(sprintf('%s Floating Wind Farm: %s ($%s$) Low/Wave Frequency Variance Balance', farm_arrangement, name, symb), 'Interpreter', 'latex')
% 
% % Loop through center turbines
% for c = 1:length(centers)
% 
%     % Get turbine and corresponding color
%     turbine = centers(c);
%     row = turbineRow(turbine, farm_arrangement);
%     colors = row_colors.(sprintf('Row%1.0f', row));
% 
%     % Plotting
%     h(c) = nexttile;
%     title(sprintf('Row %1.0f', c))
%     hold on 
% 
%     % Loop through farm spacings
%     for s = 1:length(farm_spacings)
%         farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%         available_wave_cases = fieldnames(bandPartitionedDeviations.(farm_spacing));
% 
%         % Loop through wave steepnesses
%         for st = 1:length(wave_steepnesses)
%             steep = compose('%02d', round(100 * wave_steepnesses(st)));
% 
%             % Loop through wavelengths
%             for w = 1:length(wavelengths)
%                 wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
% 
%                 % Check that wave case is available in tracking data
%                 if ismember(wave, available_wave_cases)
%                     harmonic_ratio = farm_spacings(s) / wavelengths(w);
% 
%                     % Get data
%                     total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
%                     low_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).LF^2;
%                     wave_band_variance = bandPartitionedDeviations.(farm_spacing).(wave)(turbine).(DOF).WF^2;
% 
%                     % Normalize
%                     low_norm = low_band_variance / total_variance;
%                     wave_norm = wave_band_variance / total_variance;
% 
%                     % Compute score
%                     wave_score = (wave_norm - low_norm) / (wave_norm + low_norm);
% 
%                     scatter(harmonic_ratio, wave_score, marker_size, 'filled', ...
%                             'MarkerFaceColor', colors(s,:), ...
%                             'MarkerFaceAlpha', steepness_alpha(st), ...
%                             'Marker', wave_shapes{w}, ...
%                             'HandleVisibility', 'off')
%                 end
%             end
%         end
%     end
% 
%     % Legend for each plot showing colors
%     for ww = 1:length(wavelengths)
%         plot(nan, nan, 'Color', colors(ww,:), 'linewidth', 3, ...
%             'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(ww)), 'HandleVisibility', 'on')
%     end
%     legend('Interpreter', 'latex', 'box', 'off', 'location', 'best', 'fontsize', 8)
% 
%     % Horizontal line at zero
%     yline(0, 'linestyle', '--', 'linewidth', 1, 'HandleVisibility', 'off', 'Alpha', 0.25)
%     yticks(-1:0.5:1)
%     hold off
% 
%     % Set aspect ratio equal since both values are unitless
%     % axis equal
% end
% 
% % Generate outside legend for marker shape and transparency
% axLeg = nexttile;
% hold on; axis off
% 
% % Proxies for marker shapes (spacing)
% for ww = 1:numel(wavelengths)
%     scatter(axLeg, nan, nan, marker_size, wave_shapes{ww}, 'k', 'filled', ...
%         'DisplayName', sprintf('$\\lambda = %1.0fD$', wavelengths(ww)));
% end
% 
% % Separator
% plot(axLeg, nan, nan, 'w', 'DisplayName', ' ');
% 
% % Proxies for alpha (steepness)
% for stt = 1:numel(wave_steepnesses)
%     scatter(axLeg, nan, nan, marker_size, 'o', 'k', 'filled', ...
%         'MarkerFaceAlpha', steepness_alpha(stt), ...
%         'DisplayName', sprintf('$ak = %1.2f$', wave_steepnesses(stt)));
% end
% 
% legend('Interpreter', 'latex', 'Box', 'off', 'FontSize', 10, 'Location', 'west');
% hold off
% 
% 
% % Plot axes labels
% xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
% ylabel(t, '$(P_{WF} - P_{LF}) / (P_{WF} + P_{LF})$', 'interpreter', 'latex')
% linkaxes(h, 'xy')
% xlim(h(1:length(centers)), [0.5, 2.6])
% ylim(h(1:length(centers)), [-1, 1])
% 
% clear frequency_component row frequency_name colors t d DOF name symb h s farm_spacing
% clear available_wave_cases st wave_steepness steep w wave harmonic_ratio total_variance 
% clear wave_band_variance wave_score turbine total_variance low_band_variance wave_band_variance
% clear low_norm wave_norm balance leg alpha_proxy axLeg c H legGlobal ss stt ww



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [rmsBand, psdInfo] = bandlimited_rms_fft(q, Fs, f1, f2, detrendMode)
% bandlimited_rms_fft  Compute RMS of q in the frequency band [f1,f2] using FFT
%
% This function uses Parseval's theorem to compute band-limited RMS with
% sharp frequency cutoffs (no filter roll-off losses). This ensures that
% the sum of squared RMS values across contiguous bands equals the total
% variance of the signal.
%
% Inputs:
%   q           [Nx1] signal
%   Fs          sampling frequency [Hz]
%   f1          lower band edge [Hz] (use 0 for DC)
%   f2          upper band edge [Hz] (must be <= Fs/2)
%   detrendMode 'constant' (mean removal), 'linear', or 'none'
%
% Outputs:
%   rmsBand     RMS of signal content in [f1, f2]
%   psdInfo     struct with diagnostic info (optional)

    arguments
        q (:,1) double
        Fs (1,1) double {mustBePositive}
        f1 (1,1) double {mustBeNonnegative}
        f2 (1,1) double {mustBePositive}
        detrendMode (1,:) char = 'constant'
    end

    % Input validation
    if f2 > Fs/2
        error('f2 must be <= Fs/2. Got f2=%.3f, Fs/2=%.3f', f2, Fs/2);
    end
    if f1 >= f2
        error('f1 must be < f2. Got f1=%.3f, f2=%.3f', f1, f2);
    end

    N = length(q);
    
    % Detrend signal
    switch detrendMode
        case 'constant'
            q0 = q - mean(q);
        case 'linear'
            q0 = detrend(q, 'linear');
        case 'none'
            q0 = q;
        otherwise
            q0 = q - mean(q);
    end

    % Compute FFT
    Q = fft(q0);
    
    % Frequency vector
    df = Fs / N;                    % Frequency resolution
    f = (0:N-1)' * df;              % Two-sided frequency vector
    
    % One-sided frequency vector for indexing
    if mod(N, 2) == 0
        nUnique = N/2 + 1;          % DC to Nyquist (inclusive)
    else
        nUnique = (N+1)/2;          % DC to just below Nyquist
    end
    f_onesided = f(1:nUnique);
    
    % Find frequency bins within the band [f1, f2)
    % Use >= f1 and < f2 for lower bands, include f2 if at Nyquist
    if abs(f2 - Fs/2) < df/2
        idx_band = (f_onesided >= f1) & (f_onesided <= f2);
    else
        idx_band = (f_onesided >= f1) & (f_onesided < f2);
    end
    
    % Get the magnitude squared of FFT coefficients in the band
    Q_onesided = Q(1:nUnique);
    magSq = abs(Q_onesided(idx_band)).^2;
    
    % Determine scaling for each bin
    idx_numbers = find(idx_band);
    
    power = 0;
    for i = 1:length(idx_numbers)
        k = idx_numbers(i);
        if k == 1
            % DC component (appears once)
            power = power + magSq(i);
        elseif mod(N, 2) == 0 && k == nUnique
            % Nyquist component for even N (appears once)
            power = power + magSq(i);
        else
            % All other components (appear twice: positive and negative freq)
            power = power + 2 * magSq(i);
        end
    end
    
    % Parseval's theorem: sum(|x|^2) = (1/N) * sum(|X|^2)
    variance_band = power / N^2;
    rmsBand = sqrt(variance_band);
    
    % Optional diagnostic output
    if nargout > 1
        Pxx = abs(Q_onesided).^2 / (N^2 * df);
        Pxx(2:end-1) = 2 * Pxx(2:end-1);
        if mod(N, 2) ~= 0
            Pxx(end) = 2 * Pxx(end);
        end
        
        psdInfo.f = f_onesided;
        psdInfo.Pxx = Pxx;
        psdInfo.idx = idx_band;
        psdInfo.df = df;
        psdInfo.N = N;
    end
end


