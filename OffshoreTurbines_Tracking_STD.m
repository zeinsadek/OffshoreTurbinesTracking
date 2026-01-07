%% Looking at RMS/STD of motion
% Zein Sadek

clear; close all; clc;
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")

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
% BAR CHARTS PER DOF + SINGLE STEEPNESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavelengths = [5,4,3,2];
wave_steepness = 0.12;
steep = compose('%02d', round(100 * wave_steepness));
turbine = 8;
DOF = 'y_kal';

name = names.(DOF);
symb = symbs.(DOF);

% Get proper units
if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
    units = '[Deg]';
else
    units = '[m]';
end

wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};
tmp = nan(length(farm_spacings), length(wavelengths));

% Create legend names for each wave case
legend_names = {};
for i = 1:length(wavelengths)
    wave = wavelengths(i);
    if wave ~= 0
        legend_names{i} = ['$\lambda = ', num2str(wave), 'D$'];
    else
        legend_names{i} = 'No Waves';
    end
end


% Generate array for bar chart
clc;
clear tmp
for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
    caze = strcat("WT60_", farm_spacing, "_AG0");
    fprintf('%s\n', caze)

    for w = 1:length(wavelengths)
        wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
        tmp(s,w) = deviations.(farm_spacing).(wave)(turbine).(DOF);
    end
end

% Plot
figure('color', 'white')
hb = bar(farm_spacings, tmp);
ax = gca;
box(ax, 'off');

% Color each wave case
for w = 1:length(hb)
    hb(w).FaceColor = wave_colors{w};
    hb(w).EdgeColor = 'none';
end

% Labels
% ylim([0, 10])
ylabel(sprintf('$\\sigma_{%s}$ %s', symb, units), 'interpreter', 'latex', 'fontsize', 14)
xlabel('$S_x / D$', 'interpreter', 'latex')
legend(legend_names, 'interpreter', 'latex', 'box', 'off')

title(sprintf('%s Floating Wind Farm: %s $(%s)$ Standard Deviation, $ak = %1.2f$', farm_arrangement, name, symb, wave_steepness), 'Interpreter', 'latex')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAR CHARTS PER DOF + LOOPED OVER ALL STEEPNESSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wavelengths = [5,4,3,2];
% wave_steepnesses = [0.06, 0.09, 0.12];
% turbine = 8;
% DOF = 'pitch_kal';
% wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};
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
% % Properly name DOF
% % Rotations
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
% figure('color', 'white')
% t = tiledlayout(length(wave_steepnesses), 1);
% sgtitle(sprintf('%s Floating Wind Farm: %s $(%s)$ Standard Deviation', farm_arrangement, name, symb), 'Interpreter', 'latex')
% 
% % Loop through steepnesses
% for st = 1:length(wave_steepnesses)
%     wave_steepness = wave_steepnesses(st);
%     steep = compose('%02d', round(100 * wave_steepness));
% 
%     tmp = nan(length(farm_spacings), length(wavelengths));
% 
% 
%     % Generate array for bar chart
%     clc;
%     clear tmp
%     for s = 1:length(farm_spacings)
%         farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%         caze = strcat("WT60_", farm_spacing, "_AG0");
%         fprintf('%s\n', caze)
% 
%         for w = 1:length(wavelengths)
%             wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
%             tmp(s,w) = deviations.(farm_spacing).(wave)(turbine).(DOF);
%         end
%     end
% 
%     % Plot
%     h(st) = nexttile;
%     hb = bar(farm_spacings, tmp);
%     ax = gca;
%     box(ax, 'off');
%     title(sprintf('$ak = %1.2f$', wave_steepness), 'interpreter', 'latex')
% 
%     % Color each wave case
%     for w = 1:length(hb)
%         hb(w).FaceColor = wave_colors{w};
%         hb(w).EdgeColor = 'none';
%     end
% 
%     % Legend
%     if st == 1
%         legend(legend_names, 'interpreter', 'latex', 'box', 'off', 'location', 'northoutside', 'Orientation', 'horizontal')
%     end
% 
% end
% 
% ylabel(t, sprintf('$\\sigma_{%s}$ %s', symb, units), 'interpreter', 'latex', 'fontsize', 14)
% xlabel(t, '$S_x / D$', 'interpreter', 'latex')
% linkaxes(h, 'xy')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AGAINST HARMONIC RATIO: ALL STEEPNESSES
% SINGLE DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc;
% wavelengths = [5,4,3,2];
% wave_steepnesses = [0.06, 0.09, 0.12];
% turbine = 8;
% DOF = 'yaw_kal';
% 
% % Properly name DOF
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
% clear tmp
% steepness_alpha = [0.3, 0.6, 1];
% sz = 100;
% spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
% 
% figure('color','white');
% hold on 
% for st = 1:length(wave_steepnesses)
%     wave_steepness = wave_steepnesses(st);
%     steep = compose('%02d', round(100 * wave_steepness));
%     disp(steep)
%     tmp = nan(length(farm_spacings), length(wavelengths));
%     for s = 1:length(farm_spacings)
%         farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%         caze = strcat("WT60_", farm_spacing, "_AG0");
%         fprintf('%s\n', caze)
% 
%         for w = 1:length(wavelengths)
%             wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
%             tmp(s,w) = deviations.(farm_spacing).(wave)(turbine).(DOF);
%         end
%     end
% 
% 
%     % Plot
%     for s = 1:length(farm_spacings)
%         farm_spacing_str = ['SX', num2str(farm_spacings(s) * 10)];
%         fprintf('%s\n', farm_spacing_str)
% 
%         for w = 1:length(wavelengths)
%             harmonic_ratio = farm_spacings(s) / wavelengths(w);
%             scatter(harmonic_ratio, tmp(s,w), sz, spacing_shapes{s}, 'filled', ...
%                     'MarkerFaceColor', wave_colors{w}, 'MarkerFaceAlpha', steepness_alpha(st))
%         end
%     end
% end
% 
% hold off
% ylabel(sprintf('$\\sigma_{%s}$ %s', symb, units), 'interpreter', 'latex', 'fontsize', 14)
% xlabel('$S_x / \lambda$', 'Interpreter','latex')
% xlim([0.4, 2.6])






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AGAINST HARMONIC RATIO: ALL STEEPNESSES
% LOOPED OVER ALL DOF
% NON-NORMALIZED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f', farm_arrangement, ceil(turbine / 3)), 'Interpreter', 'latex')

for d = 1:length(DOFs)
    DOF = DOFs{d};
    
    % Properly name DOF
    name = names.(DOF);
    symb = symbs.(DOF);
    
    % Get proper units
    if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
        units = '[Deg]';
    else
        units = '[m]';
    end
    

    
    %%% Plotting
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
                scatter(harmonic_ratio, deviations.(farm_spacing).(wave)(turbine).(DOF), sz, spacing_shapes{s}, 'filled', ...
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
linkaxes(h(1:3), 'y')
linkaxes(h(4:6), 'y')
xlim([0.5, 2.6])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AGAINST HARMONIC RATIO: ALL STEEPNESSES
% LOOPED OVER ALL DOF
% NORMALIZED BY AMPLITUDE AND STEEPNESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
wavelengths = [5,4,3,2];
wave_steepnesses = [0.06, 0.09, 0.12];
turbine = 8;

translations = {'x_kal', 'y_kal', 'z_kal'};
rotations = {'roll_kal', 'pitch_kal', 'yaw_kal'};

steepness_alpha = [0.3, 0.6, 1];
sz = 100;
spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};

figure('color','white')
t = tiledlayout(2,3);
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f', farm_arrangement, ceil(turbine / 3)), 'Interpreter', 'latex')

for d = 1:length(DOFs)
    DOF = DOFs{d};
    
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
    


    %%% Plotting
    clc; clear tmp
    h(d) = nexttile;
    title(sprintf('%s: $\\sigma_{%s} / %s$ [~]', name, symb, norm_symb), 'interpreter', 'latex', 'fontsize', 14)
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

                % Normalizations
                if ismember(DOF, translations)
                    amplitude = (wave_steepness * wavelengths(w) * 0.15) / (2 * pi);
                    normalization = amplitude;
                    norm_symb = 'a';
        
                elseif ismember(DOF, rotations)
                    normalization = wave_steepness;
                    norm_symb = 'ak';
                end

                % Convert rotations into rad
                if ismember(DOF, rotations)
                    scatter(harmonic_ratio, deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization, sz, spacing_shapes{s}, 'filled', ...
                            'Marker', spacing_shapes{s}, ...
                            'MarkerFaceColor', wave_colors{w}, ...
                            'MarkerFaceAlpha', steepness_alpha(st), ...
                            'HandleVisibility','off');
                % Leave translations as is
                else
                    scatter(harmonic_ratio, deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization, sz, spacing_shapes{s}, 'filled', ...
                            'Marker', spacing_shapes{s}, ...
                            'MarkerFaceColor', wave_colors{w}, ...
                            'MarkerFaceAlpha', steepness_alpha(st), ...
                            'HandleVisibility','off');
                end
            end
        end
    end



    %%% Legend (Pro shit)
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
linkaxes(h, 'xy')
xlim([0.5, 2.6])
ylim([0, 2.1])




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AGAINST HARMONIC RATIO: ALL STEEPNESSES
% LOOPED OVER ALL DOF
% NORMALIZED BY AMPLITUDE AND STEEPNESS
% LOOPED OVER ALL TURBINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc;
% wavelengths = [5,4,3,2];
% wave_steepnesses = [0.06, 0.09, 0.12];
% % turbine = 8;
% turbines_to_plot = 1:12;
% 
% % Colors per turbine
% colors.Row1 = spring(length(wavelengths));
% colors.Row2 = summer(length(wavelengths));
% colors.Row3 = autumn(length(wavelengths));
% colors.Row4 = winter(length(wavelengths));
% 
% translations = {'x_kal', 'y_kal', 'z_kal'};
% rotations = {'roll_kal', 'pitch_kal', 'yaw_kal'};
% 
% steepness_alpha = [0.3, 0.6, 1];
% sz = 100;
% spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
% wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};
% 
% figure('color','white')
% t = tiledlayout(2,3);
% sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f', farm_arrangement, ceil(turbine / 3)), 'Interpreter', 'latex')
% 
% for d = 1:length(DOFs)
%     DOF = DOFs{d};
% 
%     % Properly name DOF
%     % Rotations
%     if contains(DOF, 'pitch')
%         symb = '\theta';
%         name = 'Pitch';
%     elseif contains(DOF, 'roll')
%         symb = '\phi';
%         name = 'Roll';
%     elseif contains(DOF, 'yaw')
%         symb = '\psi';
%         name = 'Yaw';
%     % Translations
%     elseif contains(DOF, 'x_kal')
%         symb = 'x';
%         name = 'Surge';
%     elseif contains(DOF, 'y_kal')
%         symb = 'y';
%         name = 'Sway';
%     elseif contains(DOF, 'z_kal')
%         symb = 'z';
%         name = 'Heave';
%     end
% 
%     % Get proper units
%     if ismember(DOF, {'pitch_kal', 'roll_kal', 'yaw_kal'})
%         units = '[Deg]';
%     else
%         units = '[m]';
%     end
% 
%     % Normalization symboles
%     if ismember(DOF, translations)
%         norm_symb = 'a';
%     elseif ismember(DOF, rotations)
%         norm_symb = 'ak';
%     end
% 
% 
% 
%     %%% Plotting
%     clc; clear tmp
%     h(d) = nexttile;
%     title(sprintf('%s: $\\sigma_{%s} / %s$ [~]', name, symb, norm_symb), 'interpreter', 'latex', 'fontsize', 14)
%     hold on 
%     for ts = 1:length(turbines_to_plot)
%         turbine = turbines_to_plot(ts);
%         row_tag = strcat('Row', num2str(ceil(turbine / 3)));
% 
%         for st = 1:length(wave_steepnesses)
%             wave_steepness = wave_steepnesses(st);
%             steep = compose('%02d', round(100 * wave_steepness));
%             disp(steep{1})
% 
%             for s = 1:length(farm_spacings)
%                 farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%                 caze = strcat("WT60_", farm_spacing, "_AG0");
%                 fprintf('%s\n', caze)
% 
%                 for w = 1:length(wavelengths)
%                     wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
%                     harmonic_ratio = farm_spacings(s) / wavelengths(w);
% 
%                     % Normalizations
%                     if ismember(DOF, translations)
%                         amplitude = (wave_steepness * wavelengths(w) * 0.15) / (2 * pi);
%                         normalization = amplitude;
%                         norm_symb = 'a';
% 
%                     elseif ismember(DOF, rotations)
%                         normalization = wave_steepness;
%                         norm_symb = 'ak';
%                     end
% 
%                     % Convert rotations into rad
%                     if ismember(DOF, rotations)
%                         scatter(harmonic_ratio, deg2rad(deviations.(farm_spacing).(wave)(turbine).(DOF)) / normalization, sz, spacing_shapes{s}, 'filled', ...
%                                 'Marker', spacing_shapes{s}, ...
%                                 'MarkerFaceColor', colors.(row_tag)(w,:), ...
%                                 'MarkerFaceAlpha', steepness_alpha(st), ...
%                                 'HandleVisibility','off');
%                     % Leave translations as is
%                     else
%                         scatter(harmonic_ratio, deviations.(farm_spacing).(wave)(turbine).(DOF) / normalization, sz, spacing_shapes{s}, 'filled', ...
%                                 'Marker', spacing_shapes{s}, ...
%                                 'MarkerFaceColor', colors.(row_tag)(w,:), ...
%                                 'MarkerFaceAlpha', steepness_alpha(st), ...
%                                 'HandleVisibility','off');
%                     end
%                 end
%             end
%         end
%     end



    % %%% Legend (Pro shit)
    % if d == 1
    %     % Legend for color
    %     for w = 1:length(wavelengths)
    %         plot(nan, nan, 'Color', wave_colors{w}, 'linewidth', 3, ...
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
%     hold off
% end
% xlabel(t, '$S_x / \lambda$', 'Interpreter','latex')
% linkaxes(h, 'xy')
% xlim([0.5, 2.6])
% ylim([0, 2.5])












