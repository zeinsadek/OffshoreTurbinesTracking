%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFFSHORE TURBINES: STD OF DYNAMICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at RMS/STD of motion across floating wind farm
% Code plots RMS of different degrees-of-freedom and saves matfile with all
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
% COMPUTE STD OF EACH DOF: BOTH INLINE & STAGGERED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; fprintf('Computing STD...\n\n')
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
    
            % Loop through DOFs
            for d = 1:length(DOFs)
                DOF = DOFs{d};
    
                % Loop through turbines
                for t = 1:length(turbine_catalog.(arrangement).turbines)
    
                    % Compute standard deviation
                    signal = full_tracking.(arrangement).(farm_spacing).(wave)(t).(DOF);
                    deviations.(arrangement).(farm_spacing).(wave)(t).(DOF) = std(signal, 0, 'all', 'omitnan');
                end
            end
        end
    end
    fprintf('\n')
end

% Save STD to matfile for book-keeping and later coupling analysis
fprintf('\n')
save_path = fullfile(tracking_path, 'StandardDeviations', sprintf('OffshoreTracking_FullSTD_SavitskyGolay_Cutoff_%1.0f.mat', harmonic_cutoff));
fprintf('Saving Standard Deviations Data...\n')
save(save_path, 'deviations')
fprintf('Done Saving!\n\n')

% Only consider deviations from specific arrangement
deviations = deviations.(farm_arrangement);

clear s w d t a farm_spacing caze waves wave DOF save_path arrangement signal


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: NORMALIZED STD AGAINST WAVELENGTH
% LOOPED OVER ALL DOF
% CONSIDERING A SINGLE TURBINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which turbine to plot
turbine = 9;
row = turbineRow(turbine, farm_arrangement);
colors = row_colors.(sprintf('Row%1.0f', row));

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(2,3);
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f Normalized STD', farm_arrangement, row), 'Interpreter', 'latex')

% Loop through DOF
for d = 1:length(DOFs)
    DOF = DOFs{d};
    
    % Properly name DOF
    name = DOF_names.(DOF);
    symb = DOF_symbs.(DOF);
    norm_symb = getDOFnormalization(DOF);
   
    % Plotting
    h(d) = nexttile;
    title(sprintf('%s ($%s$): $\\sigma_{%s} / %s$', name, symb, symb, norm_symb), 'interpreter', 'latex', 'fontsize', 14)
    hold on 

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through wave steepesses
        % for st = 1:length(wave_steepnesses)
        for st = 1:3
            wave_steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * wave_steepness));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)
    
                    % Normalizations
                    if ismember(DOF, translations)
                        amplitude = (wave_steepness * wavelengths(w) * 0.15) / (2 * pi);
                        normalization = amplitude;
                    elseif ismember(DOF, rotations)
                        normalization = wave_steepness;
                    end
    
                    % Convert rotations into rad
                    if ismember(DOF, rotations)
                        scatter(wavelengths(w), deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization, marker_size, 'filled', ...
                                'Marker', wave_shapes{w}, ...
                                'MarkerFaceColor', colors(s,:), ...
                                'MarkerFaceAlpha', steepness_alpha(st), ...
                                'HandleVisibility','off');

                        tmpX(w) = wavelengths(w);
                        tmpY(w) = deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization;
    
                    % Leave translations as is
                    else
                        scatter(wavelengths(w), deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization, marker_size, 'filled', ...
                                'Marker', wave_shapes{w}, ...
                                'MarkerFaceColor', colors(s,:), ...
                                'MarkerFaceAlpha', steepness_alpha(st), ...
                                'HandleVisibility','off');

                        tmpX(w) = wavelengths(w);
                        tmpY(w) = deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization;
                    end
                end
            end
            
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
linkaxes(h, 'xy')
xlim([1.5, 5.5])
ylim([0, 2.2])

clear d DOF name symb norm_symb h st wave_steepness steep s caze farm_spacing tmpX tmpY P
clear w wave harmonic_ratio amplitude normalization leg row turbine t available_wave_cases


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: NORMALIZED STD AGAINST HARMONIC RATIO
% LOOPED OVER ALL DOF
% CONSIDERING A SINGLE TURBINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which turbine to plot
turbine = 9;
row = turbineRow(turbine, farm_arrangement);
colors = row_colors.(sprintf('Row%1.0f', row));

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(2,3);
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f Normalized STD', farm_arrangement, row), 'Interpreter', 'latex')

% Loop through DOF
for d = 1:length(DOFs)
    DOF = DOFs{d};
    
    % Properly name DOF
    name = DOF_names.(DOF);
    symb = DOF_symbs.(DOF);
    norm_symb = getDOFnormalization(DOF);
   
    % Plotting
    h(d) = nexttile;
    title(sprintf('%s: $\\sigma_{%s} / %s$', name, symb, norm_symb), 'interpreter', 'latex', 'fontsize', 14)
    hold on 

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through wave steepesses
        % for st = 1:length(wave_steepnesses)
        for st = 3
            wave_steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * wave_steepness));

            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)
                    harmonic_ratio = farm_spacings(s) / wavelengths(w);
    
                    % Normalizations
                    if ismember(DOF, translations)
                        amplitude = (wave_steepness * wavelengths(w) * 0.15) / (2 * pi);
                        normalization = amplitude;
                    elseif ismember(DOF, rotations)
                        normalization = wave_steepness;
                    end
    
                    % Convert rotations into rad
                    if ismember(DOF, rotations)
                        scatter(harmonic_ratio, deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization, marker_size, 'filled', ...
                                'Marker', wave_shapes{w}, ...
                                'MarkerFaceColor', colors(s,:), ...
                                'MarkerFaceAlpha', steepness_alpha(st), ...
                                'HandleVisibility','off');

                        tmpX(w) = harmonic_ratio;
                        tmpY(w) = deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization;
    
                    % Leave translations as is
                    else
                        scatter(harmonic_ratio, deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization, marker_size, 'filled', ...
                                'Marker', wave_shapes{w}, ...
                                'MarkerFaceColor', colors(s,:), ...
                                'MarkerFaceAlpha', steepness_alpha(st), ...
                                'HandleVisibility','off');

                        tmpX(w) = harmonic_ratio;
                        tmpY(w) = deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization;
                    end
                end
            end

            % Connect the dots
            plot(tmpX, tmpY, 'color', colors(s,:), 'linewidth', 1, 'HandleVisibility', 'off')
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
end
xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
linkaxes(h, 'xy')
xlim([0.5, 2.6])
ylim([0, 1.5])

clear d DOF name symb norm_symb h st wave_steepness steep s caze farm_spacing 
clear w wave harmonic_ratio amplitude normalization leg row turbine t available_wave_cases


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: NORMALIZED STD AGAINST WAVELENGTH
% LOOPED OVER ALL CENTER TURBINES
% CONSIDERING A SINGLE DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which DOF to plot
DOF = 'yaw_kal';
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);
norm_symb = getDOFnormalization(DOF);

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(1,length(centers) + 1);
sgtitle(sprintf('%s Floating Wind Farm %s ($%s$) Normalized STD: $\\sigma_{%s} / %s$', farm_arrangement, name, DOF_symbs.(DOF), symb, norm_symb), 'interpreter', 'latex', 'fontsize', 14)

% Loop through center turbines
for c = 1:length(centers)
    
    % Get turbine
    turbine = centers(c);
    colors = row_colors.(sprintf('Row%1.0f', c));
    
    % Plotting
    h(c) = nexttile;
    title(sprintf('Row %1.0f', c))
    hold on 

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through wave steepnesses
        for st = 1:length(wave_steepnesses)
            wave_steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * wave_steepness));
        
            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)

                    % Normalizations
                    if ismember(DOF, translations)
                        amplitude = (wave_steepness * wavelengths(w) * 0.15) / (2 * pi);
                        normalization = amplitude;
            
                    elseif ismember(DOF, rotations)
                        normalization = wave_steepness;
                    end
    
                    % Convert rotations into rad
                    if ismember(DOF, rotations)
                        scatter(wavelengths(w), deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization, marker_size, 'filled', ...
                                'Marker', wave_shapes{w}, ...
                                'MarkerFaceColor', colors(s,:), ...
                                'MarkerFaceAlpha', steepness_alpha(st), ...
                                'HandleVisibility','off');

                        tmpX(w) = wavelengths(w);
                        tmpY(w) = deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization;
    
                    % Leave translations as is
                    else
                        scatter(wavelengths(w), deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization, marker_size, 'filled', ...
                                'Marker', wave_shapes{w}, ...
                                'MarkerFaceColor', colors(s,:), ...
                                'MarkerFaceAlpha', steepness_alpha(st), ...
                                'HandleVisibility','off');

                        tmpX(w) = wavelengths(w);
                        tmpY(w) = deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization;
                    end
                end
            end

            % Connect the dots
            plot(tmpX, tmpY, 'color', colors(s,:), 'linewidth', 1, 'HandleVisibility', 'off')
        end
    end

    % Legend for each plot showing colors
    for ss = 1:length(farm_spacings)
        plot(nan, nan, 'Color', colors(ss,:), 'linewidth', 3, ...
            'Displayname', sprintf('$S_x = %1.1fD', farm_spacings(ss)), 'HandleVisibility', 'on')
    end
    legend('Interpreter', 'latex', 'box', 'off', 'location', 'northwest', 'fontsize', 6)
    hold off
    xticks(2:1:5)
    xlim([1.5, 5.5])
    ylim([0, 2.5])
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
linkaxes(h, 'xy')

clear d DOF name symb norm_symb h st wave_steepness steep s caze farm_spacing c axLeg ww stt
clear w wave harmonic_ratio amplitude normalization leg row turbine t available_wave_cases
clear ww 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING: NORMALIZED STD AGAINST HARMONIC RATIO
% LOOPED OVER ALL CENTER TURBINES
% CONSIDERING A SINGLE DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which DOF to plot
DOF = 'yaw_kal';
name = DOF_names.(DOF);
symb = DOF_symbs.(DOF);
norm_symb = getDOFnormalization(DOF);

% Plotting
clc; close all
figure('color','white')
t = tiledlayout(1,length(centers) + 1);
sgtitle(sprintf('%s Floating Wind Farm %s ($%s$) Normalized STD: $\\sigma_{%s} / %s$', farm_arrangement, name, DOF_symbs.(DOF), symb, norm_symb), 'interpreter', 'latex', 'fontsize', 14)

% Loop through center turbines
for c = 1:length(centers)
    
    % Get turbine
    turbine = centers(c);
    colors = row_colors.(sprintf('Row%1.0f', c));
    
    % Plotting
    h(c) = nexttile;
    title(sprintf('Row %1.0f', c))
    hold on 

    % Loop through farm spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        available_wave_cases = fieldnames(tracking.(farm_spacing));

        % Loop through wave steepnesses
        for st = 1:length(wave_steepnesses)
            wave_steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * wave_steepness));
        
            % Loop through wavelengths
            for w = 1:length(wavelengths)
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
                harmonic_ratio = farm_spacings(s) / wavelengths(w);

                % Check that wave case is available in tracking data
                if ismember(wave, available_wave_cases)

                    % Normalizations
                    if ismember(DOF, translations)
                        amplitude = (wave_steepness * wavelengths(w) * 0.15) / (2 * pi);
                        normalization = amplitude;
            
                    elseif ismember(DOF, rotations)
                        normalization = wave_steepness;
                    end
    
                    % Convert rotations into rad
                    if ismember(DOF, rotations)
                        scatter(harmonic_ratio, deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization, marker_size, 'filled', ...
                                'Marker', wave_shapes{w}, ...
                                'MarkerFaceColor', colors(s,:), ...
                                'MarkerFaceAlpha', steepness_alpha(st), ...
                                'HandleVisibility','off');
    
                    % Leave translations as is
                    else
                        scatter(harmonic_ratio, deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization, marker_size, 'filled', ...
                                'Marker', wave_shapes{w}, ...
                                'MarkerFaceColor', colors(s,:), ...
                                'MarkerFaceAlpha', steepness_alpha(st), ...
                                'HandleVisibility','off');
                    end
                end
            end
        end
    end

    % Legend for each plot showing colors
    for ss = 1:length(farm_spacings)
        plot(nan, nan, 'Color', colors(ss,:), 'linewidth', 3, ...
            'Displayname', sprintf('$S_x = %1.1fD', farm_spacings(ss)), 'HandleVisibility', 'on')
    end
    legend('Interpreter', 'latex', 'box', 'off', 'location', 'northeast', 'fontsize', 8)
    hold off
    xlim([0.5, 2.6])
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
xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
linkaxes(h, 'xy')

% ylim([0, 2])

clear d DOF name symb norm_symb h st wave_steepness steep s caze farm_spacing c axLeg ww stt
clear w wave harmonic_ratio amplitude normalization leg row turbine t available_wave_cases
clear ww 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function output = getwaveproperties(wave)
    % getwaveproperties: Reads wave case name and returns wavelength,
    % steepness, and forcing frequency

    % Split wave name at underscore
    tmp = split(wave, '_');
    wavelength_string = tmp{1};
    steepness_string = tmp{2};

    % Get wavelengths and steepness from strings
    wavelength = wavelength_from_string.(wavelength_string);
    steepness = steepness_from_string.(steepness_string);

    % Get wave amplitude from steepness and wavelength
    amplitude = (steepness * wavelength * 0.15) / (2 * pi);

    % Returned values
    output.wavelength = wavelength;
    output.steepness = steepness;
    output.frequency = forcing_frequencies.(wavelength);
    output.amplitude = amplitude;
end



function unit = getDOFunits(DOF)
    % getDOFunits: scans input DOF string and outputs the corresponding
    % units

    % Get proper units
    if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
        unit = '[Deg]';
    elseif ismember(DOF, {'x_kal', 'y_kal', 'x_kal'})
        unit = '[m]';
    end
end



function norm_symb = getDOFnormalization(DOF)
    % getDOFnormalization: scans input DOF string and outputs the
    % corresponding normalizatio unit

    if ismember(DOF, {'x_kal', 'y_kal', 'z_kal'})
        norm_symb = 'a';
    elseif ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
        norm_symb = 'ak';
    end
end




