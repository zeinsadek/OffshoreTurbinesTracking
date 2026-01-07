%% Plotting signals for graphics
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
% PLOT SIGNALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

farm_spacing = 'SX50';
wave = 'LM5_AK12';
turbine = 8;

colors = cool(length(DOFs));

% Motion names
names.x_kal = 'Surge';
names.y_kal = 'Heave';
names.z_kal = 'Sway';
names.roll_kal = 'Roll';
names.pitch_kal = 'Pitch';
names.yaw_kal = 'Yaw';

% Motion symbols
symbs.x_kal = 'x';
symbs.y_kal = 'y';
symbs.z_kal = 'z';
symbs.roll_kal = '\phi';
symbs.pitch_kal = '\theta';
symbs.yaw_kal = '\psi';

figure('color', 'white')
t = tiledlayout(length(DOFs), 1);

for d = 1:length(DOFs)
    DOF = DOFs{d};
    h(d) = nexttile;
    plot(tracking.(farm_spacing).(wave)(turbine).time, tracking.(farm_spacing).(wave)(turbine).(DOF), ...
         'linewidth', 1, 'color', colors(d,:))
    set(gca, 'box', 'off')

    % Make it sexy
    if d ~= length(DOFs)
        ax = gca;
        ax.XAxis.Visible = 'off';
        set(ax,'xtick',[])
    end

    % Axes limits
    if ismember(d, [1,2,3])
        ylim([-0.02, 0.02])
    else
        ylim([-30, 30])
    end


    
    title(sprintf('%s $\\left( %s \\right)$', names.(DOF), symbs.(DOF)), 'interpreter', 'latex')

end

linkaxes(h(1:3), 'y')
linkaxes(h(4:6), 'y')
xlabel(t, 'Time [s]')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT STD BAR CHART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

farm_spacing = 'SX50';
wave = 'LM5_AK12';
turbine = 8;

colors = cool(length(DOFs));

% Motion names
names.x_kal = 'Surge';
names.y_kal = 'Heave';
names.z_kal = 'Sway';
names.roll_kal = 'Roll';
names.pitch_kal = 'Pitch';
names.yaw_kal = 'Yaw';

% Motion symbols
symbs.x_kal = 'x';
symbs.y_kal = 'y';
symbs.z_kal = 'z';
symbs.roll_kal = '\phi';
symbs.pitch_kal = '\theta';
symbs.yaw_kal = '\psi';

clear tmp labels
for d = 1:length(DOFs)
    DOF = DOFs{d};
    tmp(d) = deviations.(farm_spacing).(wave)(turbine).(DOF);
    labels{d} = sprintf('%s $\\left( %s \\right)$', names.(DOF), symbs.(DOF));

end

figure('color', 'white')
tiledlayout(1,2)

% Translations
nexttile
b = bar(labels(1:3), tmp(1:3));
b.FaceColor = 'flat';
b.CData     = colors(1:3,:);   % one RGB row per bar
set(gca,"TickLabelInterpreter",'latex')
ylabel('$\sigma$ [m]', 'interpreter', 'latex')

% Rotations
nexttile
b = bar(labels(4:6), tmp(4:6));
b.FaceColor = 'flat';
b.CData     = colors(4:6,:);
set(gca,"TickLabelInterpreter",'latex')
ylabel('$\sigma$ [Deg]', 'interpreter', 'latex')







