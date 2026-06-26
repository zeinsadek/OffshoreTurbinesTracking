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
DOF_symbs.x_kal = '$\Delta x$';
DOF_symbs.y_kal = '$\Delta y$';
DOF_symbs.z_kal = '$\Delta z$';
DOF_symbs.roll_kal = '$\phi$';
DOF_symbs.pitch_kal = '$\theta$';
DOF_symbs.yaw_kal = '$\psi$';
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



%% Tes with all turbines (fluctuations)

tickFontSize = 8;
labelFontSize = 10;

linewidth = 0.5;
colors = slanCM('dusk', length(DOFs));


farm_arrangement = 'Inline';
farm_spacing = 'SX50';
wave = 'LM5_AK12';
turbine = 2;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 7, 5]);
tiledlayout(length(DOFs), 1, 'padding', 'none', 'TileSpacing', 'none')
for i = 1:length(DOFs)
    DOF = DOFs{i};
    h(i) = nexttile;

    if i <= 3
        scale = 1E3;
    else
        scale = 1;
    end

    hold on
    plot(tracking.(farm_spacing).(wave)(turbine).time, ...
         tracking.(farm_spacing).(wave)(turbine).(DOF) * scale, ...
         'linewidth', linewidth, 'color', colors(i,:))
    hold off

    set(h(i), 'FontSize', tickFontSize)
    set(h(i), 'TickLabelInterpreter', 'latex');

    if i > 3
         label = sprintf('%s\\,\\,\\,', DOF_symbs.(DOF));
    else
         label = sprintf('%s\\,', DOF_symbs.(DOF));
    end

    % if i == 5
    %     label = sprintf('%s\\,\\,\\,', DOF_symbs.(DOF));
    % end

    ylabel(label, 'FontSize', labelFontSize, 'interpreter', 'latex', 'rotation', 0)
    yticks([])

    if i < 6
        ax = gca;
        ax.XAxis.Visible = 'off';
    end
end

linkaxes(h, 'xy')
xlim([0, 60])

xlabel('Time [s]', 'interpreter', 'latex', 'fontsize', labelFontSize)


% Test save
save_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/Publishing/PowerMotion/Tracking';
save_name = 'Offshore_MotionPower_TrackingSampleSignals.pdf';
pause(5)
exportgraphics(fig, fullfile(save_folder, save_name), 'resolution', 600)
close all