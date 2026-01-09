%% Looking at band-limited RMS/STD of motion (CORRECTED VERSION)
% Zein Sadek
% 
% Key changes from original:
% 1. Uses FFT-based band-limited RMS (no filter roll-off losses)
% 2. Bands are defined to be perfectly contiguous: [0, 0.5*fw], [0.5*fw,
% 1.5*fw], [1.5*fw, Fs/2]
% 3. Includes verification that sum of squared band RMS equals total variance

clear; close all; clc;
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT TRACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

offshore_path = "/Users/zeinsadek/Desktop/Experiments/Offshore";
power_path = fullfile(offshore_path, "Power/Data/Matfiles");
tracking_path = fullfile(offshore_path, "Tracking/Data/Matfiles");

farm_arrangement = "Staggered";

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STD OF EACH DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DOFs = {'x_kal', 'y_kal', 'z_kal', 'roll_kal', 'pitch_kal', 'yaw_kal'};

clc;
for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
    caze = strcat("WT60_", farm_spacing, "_AG0");
    fprintf('%s\n', caze)
    
    tmp = tracking.(farm_spacing);
    waves = fieldnames(tmp);

    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};

        % Loop through DOFs
        for d = 1:length(DOFs)
            DOF = DOFs{d};

            % Loop through turbines
            for t = 1:length(turbines)
                % Compute standard deviation
                deviations.(farm_spacing).(wave)(t).(DOF) = std(tracking.(farm_spacing).(wave)(t).(DOF), 0, 'all', 'omitnan');
            end
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE BAND-LIMITED RMS OF EACH DOF (FFT METHOD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% KEY CHANGES:
% 1. Using FFT-based method for sharp frequency cutoffs
% 2. Bands are perfectly contiguous: [0, 0.5*fw], [0.5*fw, 1.5*fw], [1.5*fw, Fs/2]
% 3. This ensures sum of squared RMS values = total variance

Fs = 30;
DOFs = {'x_kal', 'y_kal', 'z_kal', 'roll_kal', 'pitch_kal', 'yaw_kal'};

clc;
fprintf('Computing band-limited RMS using FFT method...\n\n');

for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
    caze = strcat("WT60_", farm_spacing, "_AG0");
    fprintf('%s\n', caze)
    
    tmp = tracking.(farm_spacing);
    waves = fieldnames(tmp);

    % Loop through waves
    for w = 1:length(waves)
        % Get wave and wave frequency
        wave = waves{w};
        split_wave = split(wave, '_');
        wavelength_name = split_wave{1};

        if ~strcmp(wavelength_name, 'LM0')
            fw = forcing_frequencies.(wavelength_name);

            % Define contiguous bands spanning [0, Fs/2]
            % No gaps, no overlaps
            bands.LF   = [0,       0.5*fw];     % DC to half wave frequency
            bands.Wave = [0.5*fw,  1.5*fw];     % Wave frequency band
            bands.HF   = [1.5*fw,  Fs/2];       % High frequency to Nyquist
    
            % Loop through DOFs
            for d = 1:length(DOFs)
                DOF = DOFs{d};
    
                % Loop through turbines
                for t = 1:length(turbines)
                    % Signal
                    q = tracking.(farm_spacing).(wave)(t).(DOF);
    
                    % Compute low frequency portion (includes DC)
                    rLF = bandlimited_rms_fft(q, Fs, bands.LF(1), bands.LF(2), 'constant');
                    
                    % Compute wave frequency portion
                    rWave = bandlimited_rms_fft(q, Fs, bands.Wave(1), bands.Wave(2), 'constant');
    
                    % Compute high frequency portion (up to Nyquist)
                    rHF = bandlimited_rms_fft(q, Fs, bands.HF(1), bands.HF(2), 'constant');

                    % Save
                    bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).LF = rLF;
                    bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).WF = rWave;
                    bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).HF = rHF;
    
                end
            end
        end
    end
end

fprintf('\nDone!\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERIFY: Sum of squared bands should equal total variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
fprintf('=== VERIFICATION: Parseval''s Theorem Check ===\n\n');

spacing = 'SX50';
wave = 'LM5_AK12';
turbine = 1;
DOF = 'yaw_kal';

full = deviations.(spacing).(wave)(turbine).(DOF);
LF = bandfilteredDeviations.(spacing).(wave)(turbine).(DOF).LF;
WF = bandfilteredDeviations.(spacing).(wave)(turbine).(DOF).WF;
HF = bandfilteredDeviations.(spacing).(wave)(turbine).(DOF).HF;

% Sum of squares of partitioned signals
ss = LF^2 + WF^2 + HF^2;
total_squared = full^2;

fprintf('Case: %s, %s, Turbine %d, DOF: %s\n', spacing, wave, turbine, DOF);
fprintf('--------------------------------------------\n');
fprintf('Band RMS values:\n');
fprintf('  LF:   %.6f\n', LF);
fprintf('  Wave: %.6f\n', WF);
fprintf('  HF:   %.6f\n', HF);
fprintf('\nFull signal STD: %.6f\n', full);
fprintf('\nSum of squared bands:  %.6f\n', ss);
fprintf('Squared full STD:      %.6f\n', total_squared);
fprintf('\nRatio (should be ~1.0): %.6f\n', ss / total_squared);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK ALL CASES - Find worst closure errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
fprintf('=== Checking closure errors for all cases ===\n\n');

ratios = [];
worst_cases = struct('ratio', [], 'spacing', {}, 'wave', {}, 'turbine', [], 'DOF', {});

for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
    tmp = tracking.(farm_spacing);
    waves = fieldnames(tmp);

    for w = 1:length(waves)
        wave = waves{w};
        split_wave = split(wave, '_');
        wavelength_name = split_wave{1};

        if ~strcmp(wavelength_name, 'LM0') && isfield(bandfilteredDeviations.(farm_spacing), wave)
            for d = 1:length(DOFs)
                DOF = DOFs{d};
                for t = 1:length(turbines)
                    full = deviations.(farm_spacing).(wave)(t).(DOF);
                    LF = bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).LF;
                    WF = bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).WF;
                    HF = bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).HF;
                    
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

fprintf('Closure ratio statistics:\n');
fprintf('  Min:    %.6f\n', min(ratios));
fprintf('  Max:    %.6f\n', max(ratios));
fprintf('  Mean:   %.6f\n', mean(ratios));
fprintf('  Median: %.6f\n', median(ratios));
fprintf('  Std:    %.6f\n', std(ratios));

fprintf('\nCases with >1%% error: %d out of %d\n', length(worst_cases), length(ratios));

if ~isempty(worst_cases)
    fprintf('\nWorst cases:\n');
    [~, sortIdx] = sort(abs([worst_cases.ratio] - 1), 'descend');
    for i = 1:min(5, length(sortIdx))
        wc = worst_cases(sortIdx(i));
        fprintf('  %s, %s, T%d, %s: ratio = %.4f\n', ...
            wc.spacing, wc.wave, wc.turbine, wc.DOF, wc.ratio);
    end
end

% Histogram
% figure;
% histogram(ratios, 50);
% xlabel('Closure Ratio (sum of squared bands / total variance)');
% ylabel('Count');
% title('Distribution of Parseval''s Theorem Closure Ratios');
% xline(1, 'r--', 'LineWidth', 2);
% legend('Ratios', 'Perfect closure');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOMENCLATURE FOR DOFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DOF names
names.x_kal = 'Surge';
names.y_kal = 'Heave';
names.z_kal = 'Sway';
names.roll_kal = 'Roll';
names.pitch_kal = 'Pitch';
names.yaw_kal = 'Yaw';

% DOF symbols
symbs.x_kal = 'x';
symbs.y_kal = 'y';
symbs.z_kal = 'z';
symbs.roll_kal = '\phi';
symbs.pitch_kal = '\theta';
symbs.yaw_kal = '\psi';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: STACKED BAR CHART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequencies = {'LF', 'WF', 'HF'};
% wavelengths = [5,4,3,2];
% wave_steepness = 0.12;
% steep = compose('%02d', round(100 * wave_steepness));
% turbine = 8;
% DOF = 'z_kal';
% 
% name = names.(DOF);
% symb = symbs.(DOF);
% 
% % Get proper units
% if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
%     units = '[Deg]';
% else
%     units = '[m]';
% end
% 
% wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};
% tmp = nan(length(farm_spacings), length(wavelengths), 3);
% 
% % Create legend names for each wave case
% legend_names = {};
% for i = 1:length(wavelengths)
%     wave = wavelengths(i);
%     if wave ~= 0
%         legend_names{i} = ['$\lambda = ', num2str(wave), 'D$'];
%     else
%         legend_names{i} = 'No Waves';
%     end
% end
% 
% 
% % Generate array for bar chart
% clc;
% for s = 1:length(farm_spacings)
%     farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%     caze = strcat("WT60_", farm_spacing, "_AG0");
%     fprintf('%s\n', caze)
% 
%     for w = 1:length(wavelengths)
%         wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
% 
%         for f = 1:3
%             tmp(s,w,f) = bandfilteredDeviations.(farm_spacing).(wave)(turbine).(DOF).(frequencies{f});
%         end
%     end
% end
% 
% 
% stackAlpha  = linspace(0.4,1,size(tmp,3));  % opacity per STACK
% groupLabels = arrayfun(@(x) sprintf("%g", x), farm_spacings, 'UniformOutput', false);
% 
% 
% % Sort by xVals (or use your own custom order)
% [xSorted, idx] = sort(farm_spacings, 'ascend');
% 
% tmp  = tmp(idx,:,:);
% labels = groupLabels(idx);
% 
% 
% figure('color', 'white')
% hold on
% plotBarStackGroups(tmp, wave_colors, stackAlpha, labels);
% 
% ylabel(sprintf('$\\sigma_{%s}$ %s', symb, units), 'interpreter', 'latex', 'fontsize', 14)
% xlabel('$S_x / D$', 'interpreter', 'latex')
% legend(legend_names, 'interpreter', 'latex', 'box', 'off')
% 
% title(sprintf('%s Floating Wind Farm: %s $(%s)$ Standard Deviation, $ak = %1.2f$', farm_arrangement, name, symb, wave_steepness), 'Interpreter', 'latex')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AGAINST HARMONIC RATIO: ALL STEEPNESSES
% LOOPED OVER ALL DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the variance associated with the wave-frequency normalized
% by the total variance
% Represents the percentage of the energy contained in this band, and varys
% from 0 to 1

clc;
wavelengths = [5,4,3,2];
wave_steepnesses = [0.06, 0.09, 0.12];
turbine = 11;

freq = 'HF';

% Titles based on frequency
if strcmp(freq, 'HF')
    freq_name = 'High Frequency';
elseif strcmp(freq, 'WF')
    freq_name = 'Wave Frequency';
elseif strcmp(freq, 'LF')
    freq_name = 'Low Frequency';
end


steepness_alpha = [0.3, 0.6, 1];
sz = 100;
spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};

clc; close all
figure('color','white')
t = tiledlayout(2,3);
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f %s Wave Score', farm_arrangement, ceil(turbine / 3), freq_name), 'Interpreter', 'latex')

for d = 1:length(DOFs)
    DOF = DOFs{d};
    name = names.(DOF);
    symb = symbs.(DOF);
    
    % Get proper units
    if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
        units = '[Deg]';
    else
        units = '[m]';
    end
    

    % Plotting
    clc; clear tmp
    h(d) = nexttile;
    title(sprintf('%s: $\\sigma_{%s}$ %s', name, symb, units), 'interpreter', 'latex', 'fontsize', 14)
    hold on 
    for st = 1:length(wave_steepnesses)
        wave_steepness = wave_steepnesses(st);
        steep = compose('%02d', round(100 * wave_steepness));
        disp(steep{1})
        
        for s = 1:length(farm_spacings)
            farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
            caze = strcat("WT60_", farm_spacing, "_AG0");
            fprintf('%s\n', caze)
        
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
                harmonic_ratio = farm_spacings(s) / wavelengths(w);

                total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
                wave_band_variance = bandfilteredDeviations.(farm_spacing).(wave)(turbine).(DOF).(freq)^2;
                wave_score = wave_band_variance ./ total_variance;

                scatter(harmonic_ratio, wave_score, sz, spacing_shapes{s}, 'filled', ...
                        'MarkerFaceColor', wave_colors{w}, 'MarkerFaceAlpha', steepness_alpha(st), ...
                        'HandleVisibility', 'off')
            end
        end
    end



    %%% Legend
    if d == 1
        % Legend for color
        for w = 1:length(wavelengths)
            plot(nan, nan, 'Color', wave_colors{w}, 'linewidth', 3, ...
                'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for s = 1:length(farm_spacings)
            scatter(nan, nan, sz, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker alpha
        for st = 1:length(wave_steepnesses)
            scatter(nan, nan, sz, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
                    'markerfacealpha', steepness_alpha(st), ...
                    'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
        end

        leg = legend('interpreter', 'latex', 'box', 'off');
        leg.Layout.Tile = 'east';
    end
    hold off
end
xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
ylabel(t, sprintf('$\\sigma_{%s}^2 / \\sigma^2$', freq), 'interpreter', 'latex')
linkaxes(h, 'xy')
xlim([0.5, 2.6])
ylim([0, 1])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AGAINST HARMONIC RATIO: ALL STEEPNESSES
% LOOPED OVER ALL DOF
% SINGLE TURBINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the balance between the low and wave frequency contribution

clc;
wavelengths = [5,4,3,2];
wave_steepnesses = [0.06, 0.09, 0.12];
turbine = 8;

steepness_alpha = [0.3, 0.6, 1];
sz = 100;
spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};

figure('color','white')
t = tiledlayout(2,3);
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f Low/Wave Frequency Balance', farm_arrangement, ceil(turbine / 3)), 'Interpreter', 'latex')

for d = 1:length(DOFs)
    DOF = DOFs{d};
    name = names.(DOF);
    symb = symbs.(DOF);
   
    % Get proper units
    if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
        units = '[Deg]';
    else
        units = '[m]';
    end
    
    
    % Plotting
    clc; clear tmp
    h(d) = nexttile;
    title(sprintf('%s: $\\sigma_{%s}$ %s', name, symb, units), 'interpreter', 'latex', 'fontsize', 14)
    hold on 
    for st = 1:length(wave_steepnesses)
        wave_steepness = wave_steepnesses(st);
        steep = compose('%02d', round(100 * wave_steepness));
        disp(steep{1})
        
        for s = 1:length(farm_spacings)
            farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
            caze = strcat("WT60_", farm_spacing, "_AG0");
            fprintf('%s\n', caze)
        
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
                harmonic_ratio = farm_spacings(s) / wavelengths(w);

                % Get data
                total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
                low_band_variance = bandfilteredDeviations.(farm_spacing).(wave)(turbine).(DOF).LF^2;
                wave_band_variance = bandfilteredDeviations.(farm_spacing).(wave)(turbine).(DOF).WF^2;

                % Normalize
                low_norm = low_band_variance / total_variance;
                wave_norm = wave_band_variance / total_variance;

                % Compute score
                wave_score = (wave_norm - low_norm) / (wave_norm + low_norm);

                scatter(harmonic_ratio, wave_score, sz, spacing_shapes{s}, 'filled', ...
                        'MarkerFaceColor', wave_colors{w}, 'MarkerFaceAlpha', steepness_alpha(st), ...
                        'HandleVisibility', 'off')
            end
        end
    end



    %%% Legend
    if d == 1
        % Legend for color
        for w = 1:length(wavelengths)
            plot(nan, nan, 'Color', wave_colors{w}, 'linewidth', 3, ...
                'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker shape
        for s = 1:length(farm_spacings)
            scatter(nan, nan, sz, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
        end

        % White space
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
        plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')

        % Legend for marker alpha
        for st = 1:length(wave_steepnesses)
            scatter(nan, nan, sz, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
                    'markerfacealpha', steepness_alpha(st), ...
                    'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
        end

        leg = legend('interpreter', 'latex', 'box', 'off');
        leg.Layout.Tile = 'east';
    end
    hold off
end
xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
ylabel(t, '$P_{WF} - P_{LF} / P_{WF} + P_{LF}$', 'interpreter', 'latex')
linkaxes(h, 'xy')
xlim([0.5, 2.6])
ylim([-1, 1])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AGAINST HARMONIC RATIO: ALL STEEPNESSES
% LOOPED OVER CENTER TURBINES
% SINGLE DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the balance between the low and wave frequency contribution

% SX30, LM5_AK12 case has bad data for staggered. 
% Consider removing from plots


clc; close all
wavelengths = [5,4,3,2];
wave_steepnesses = [0.06, 0.09, 0.12];

if strcmp(farm_arrangement, 'Inline')
    centers = [2,5,8,10];
else
    centers = [2,4,7,9];
end

DOF = 'yaw_kal';

translations = {'x_kal', 'y_kal', 'z_kal'};
rotations = {'roll_kal', 'pitch_kal', 'yaw_kal'};

% Properly name DOF
name = names.(DOF);
symb = symbs.(DOF);

% Get proper units
if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
    units = '[Deg]';
else
    units = '[m]';
end

% Normalization symboles
if ismember(DOF, translations)
    norm_symb = 'a';
elseif ismember(DOF, rotations)
    norm_symb = 'ak';
end

% Colors per row
row_colors.Row1 = flipud(slanCM(45, 2 * length(wavelengths)));
row_colors.Row2 = flipud(slanCM(47, 2 * length(wavelengths)));
row_colors.Row3 = flipud(slanCM(34, 2 * length(wavelengths)));
row_colors.Row4 = flipud(slanCM(35, 2 * length(wavelengths)));

steepness_alpha = [0.3, 0.6, 1];
sz = 100;
spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};

figure('color','white')
t = tiledlayout(1,length(centers));
sgtitle(sprintf('%s Floating Wind Farm: %s ($\\sigma_{%s}$)  Low/Wave Frequency Balance', farm_arrangement, name, symb), 'Interpreter', 'latex')

for c = 1:length(centers)

    turbine = centers(c);
    colors = row_colors.(sprintf('Row%1.0f', c));

    % Plotting
    clc; clear tmp
    h(c) = nexttile;
    % title(sprintf('%s: $\\sigma_{%s}$ %s', name, symb, units), 'interpreter', 'latex', 'fontsize', 14)
    title(sprintf('Row%1.0f', c))
    hold on 
    for st = 1:length(wave_steepnesses)
        wave_steepness = wave_steepnesses(st);
        steep = compose('%02d', round(100 * wave_steepness));
        disp(steep{1})
        
        for s = 1:length(farm_spacings)
            farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
            caze = strcat("WT60_", farm_spacing, "_AG0");
            fprintf('%s\n', caze)
        
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
                harmonic_ratio = farm_spacings(s) / wavelengths(w);

                % Get data
                total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
                low_band_variance = bandfilteredDeviations.(farm_spacing).(wave)(turbine).(DOF).LF^2;
                wave_band_variance = bandfilteredDeviations.(farm_spacing).(wave)(turbine).(DOF).WF^2;

                % Normalize
                low_norm = low_band_variance / total_variance;
                wave_norm = wave_band_variance / total_variance;

                % Compute score
                wave_score = (wave_norm - low_norm) / (wave_norm + low_norm);

                scatter(harmonic_ratio, wave_score, sz, spacing_shapes{s}, 'filled', ...
                        'MarkerFaceColor', colors(w,:), 'MarkerFaceAlpha', steepness_alpha(st), ...
                        'HandleVisibility', 'off')
            end
        end
    end



    % %%% Legend
    % if d == 1
    %     % Legend for color
    %     for w = 1:length(wavelengths)
    %         plot(nan, nan, 'Color', colors(w,:), 'linewidth', 3, ...
    %             'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
    %     end
    % 
    %     % White space
    %     plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
    %     plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
    % 
    %     % Legend for marker shape
    %     for s = 1:length(farm_spacings)
    %         scatter(nan, nan, sz, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
    %                 'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
    %     end
    % 
    %     % White space
    %     plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
    %     plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
    % 
    %     % Legend for marker alpha
    %     for st = 1:length(wave_steepnesses)
    %         scatter(nan, nan, sz, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
    %                 'markerfacealpha', steepness_alpha(st), ...
    %                 'Displayname', sprintf('$ak = %1.2f$', wave_steepnesses(st)))
    %     end
    % 
    %     leg = legend('interpreter', 'latex', 'box', 'off');
    %     leg.Layout.Tile = 'east';
    % end
    yticks(-1:0.5:1)
    hold off
end
xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
ylabel(t, '$P_{WF} - P_{LF} / P_{WF} + P_{LF}$', 'interpreter', 'latex')
linkaxes(h, 'xy')
xlim([0.5, 2.6])
ylim([-1, 1])





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AGAINST HARMONIC RATIO: ALL STEEPNESSES
% LOOPED OVER ALL DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting here the balance between the low and high frequency contribution

% clc;
% wavelengths = [5,4,3,2];
% wave_steepnesses = [0.06, 0.09, 0.12];
% turbine = 8;
% 
% steepness_alpha = [0.3, 0.6, 1];
% sz = 100;
% spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
% wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};
% 
% figure('color','white')
% t = tiledlayout(2,3);
% sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f Low/High Frequency Balance', farm_arrangement, ceil(turbine / 3)), 'Interpreter', 'latex')
% 
% for d = 1:length(DOFs)
%     DOF = DOFs{d};
%     name = names.(DOF);
%     symb = symbs.(DOF);
% 
%     % Get proper units
%     if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
%         units = '[Deg]';
%     else
%         units = '[m]';
%     end
% 
% 
%     % Plotting
%     clc; clear tmp
%     h(d) = nexttile;
%     title(sprintf('%s: $\\sigma_{%s}$ %s', name, symb, units), 'interpreter', 'latex', 'fontsize', 14)
%     hold on 
%     for st = 1:length(wave_steepnesses)
%         wave_steepness = wave_steepnesses(st);
%         steep = compose('%02d', round(100 * wave_steepness));
%         disp(steep{1})
% 
%         for s = 1:length(farm_spacings)
%             farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%             caze = strcat("WT60_", farm_spacing, "_AG0");
%             fprintf('%s\n', caze)
% 
%             for w = 1:length(wavelengths)
%                 wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
%                 harmonic_ratio = farm_spacings(s) / wavelengths(w);
% 
%                 % Get data
%                 total_variance = deviations.(farm_spacing).(wave)(turbine).(DOF)^2;
%                 low_band_variance = bandfilteredDeviations.(farm_spacing).(wave)(turbine).(DOF).LF^2;
%                 high_band_variance = bandfilteredDeviations.(farm_spacing).(wave)(turbine).(DOF).HF^2;
% 
%                 % Normalize
%                 low_norm = low_band_variance / total_variance;
%                 high_norm = high_band_variance / total_variance;
% 
%                 % Compute score
%                 wave_score = (high_norm - low_norm) / (high_norm + low_norm);
% 
%                 scatter(harmonic_ratio, wave_score, sz, spacing_shapes{s}, 'filled', ...
%                         'MarkerFaceColor', wave_colors{w}, 'MarkerFaceAlpha', steepness_alpha(st), ...
%                         'HandleVisibility', 'off')
%             end
%         end
%     end
% 
% 
% 
%     %%% Legend
%     if d == 1
%         % Legend for color
%         for w = 1:length(wavelengths)
%             plot(nan, nan, 'Color', wave_colors{w}, 'linewidth', 3, ...
%                 'Displayname', sprintf('$\\lambda = %1.0fD$', wavelengths(w)), 'HandleVisibility', 'on')
%         end
% 
%         % White space
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
% 
%         % Legend for marker shape
%         for s = 1:length(farm_spacings)
%             scatter(nan, nan, sz, spacing_shapes{s}, 'black', 'filled', 'HandleVisibility', 'on', ...
%                     'DisplayName', sprintf('$S_x = %1.1fD', farm_spacings(s)))
%         end
% 
%         % White space
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
%         plot(nan, nan, 'color', 'white', 'HandleVisibility', 'on', 'displayname', '')
% 
%         % Legend for marker alpha
%         for st = 1:length(wave_steepnesses)
%             scatter(nan, nan, sz, 'o', 'black', 'filled', 'HandleVisibility', 'on', ...
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
% ylabel(t, '$P_{HF} - P_{LF} / P_{HF} + P_{LF}$', 'interpreter', 'latex')
% linkaxes(h, 'xy')
% xlim([0.5, 2.6])
% ylim([-1, 1])
















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT-BASED BAND-LIMITED RMS FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



function [] = plotBarStackGroups(stackData, barColors, stackAlpha, groupLabels)
% Plot a set of stacked bars, but group them according to labels provided.
%
% Params: 
%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)

NumGroupsPerAxis = size(stackData, 1);
NumStacksPerGroup = size(stackData, 2);


% Count off the number of bins
groupBins = 1:NumGroupsPerAxis;
MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
groupOffset = MaxGroupWidth/NumStacksPerGroup;
figure('color', 'white')
hold on
for i=1:NumStacksPerGroup

    Y = squeeze(stackData(:,i,:));
    
    % Center the bars:
    internalPosCount = i - ((NumStacksPerGroup+1) / 2);
    
    % Offset the group draw positions:
    groupDrawPos = (internalPosCount)* groupOffset + groupBins;
    
    h(i,:) = bar(Y, 'stacked', 'HandleVisibility', 'off');
    set(h(i,:),'BarWidth',groupOffset);
    set(h(i,:),'XData',groupDrawPos);
    
    % ---- color stacks within THIS bar ----
    for k = 1:size(Y,2)   % stack elements
        % h(i,k).FaceColor = barColors(i,:);
        h(i,k).FaceColor = barColors{i};
        h(i,k).FaceAlpha = stackAlpha(k);
        h(i,k).EdgeColor = 'none';
    end

end
hold off
set(gca,'XTickMode','manual');
set(gca,'XTick',1:NumGroupsPerAxis);
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',groupLabels);
end 