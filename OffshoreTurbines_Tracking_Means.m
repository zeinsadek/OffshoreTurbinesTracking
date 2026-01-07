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
% COMPUTE MEAN OF EACH DOF
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
                averages.(farm_spacing).(wave)(t).(DOF) = mean(tracking.(farm_spacing).(wave)(t).(DOF), 'all', 'omitnan');
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
% PLOTTING MEAN OF ALL DOF AGAINST HARMONIC RATIO
% SINGLE TURBINE
% PLOTTED AS A TILEDLAYOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
wavelengths = [5,4,3,2];
wave_steepnesses = [0.06, 0.09, 0.12];
turbine = 8;

steepness_alpha = [0.3, 0.6, 1];
sz = 100;
spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};

clc; close all
figure('color','white')
t = tiledlayout(2,3);
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f Mean Dynamics', farm_arrangement, ceil(turbine / 3)), 'Interpreter', 'latex')

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
                average = averages.(farm_spacing).(wave)(turbine).(DOF);

                scatter(harmonic_ratio, average, sz, spacing_shapes{s}, 'filled', ...
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
linkaxes(h, 'xy')
xlim([0.5, 2.6])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING MEAN OF PITCH AGAINST HARMONIC RATIO
% ALL CENTER TURBINES
% PLOTTED AS A TILEDLAYOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
wavelengths = [5,4,3,2];
wave_steepnesses = [0.06, 0.09, 0.12];
centers = [2,5,8,11];
DOF = 'pitch_kal';

steepness_alpha = [0.3, 0.6, 1];
sz = 100;
spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};

clc; close all
figure('color','white')
t = tiledlayout(1,length(centers));
sgtitle(sprintf('%s Floating Wind Farm: Row %1.0f Mean Dynamics', farm_arrangement, ceil(turbine / 3)), 'Interpreter', 'latex')

for c = 1:length(centers)

    turbine = centers(c);
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
    h(c) = nexttile;
    title(sprintf('Row %1.0f %s', c, name), 'interpreter', 'latex', 'fontsize', 14)
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
                average = averages.(farm_spacing).(wave)(turbine).(DOF);

                scatter(harmonic_ratio, average, sz, spacing_shapes{s}, 'filled', ...
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
linkaxes(h, 'xy')
xlim([0.5, 2.6])


