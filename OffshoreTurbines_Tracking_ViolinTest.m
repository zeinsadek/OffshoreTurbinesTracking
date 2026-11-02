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



%% Tes with all turbines (fluctuations)


% full_tracking.(arrangement).(farm_spacing).(wave)(t).(DOF)

DOF = 'pitch_kal';
farm_arrangement = 'Inline';
farm_spacing = 'SX50';
wave = 'LM5_AK12';

clc; close all
figure('color', 'white')
hold on
for i = 1:3:12
    tmp = full_tracking.(farm_arrangement).(farm_spacing).(wave)(i).(DOF);
    tmp = tmp - mean(tmp, 'all', 'omitnan');
    violinplot(i, tmp)
end
hold off
yline(0)
xticks(1:12)
xlabel('Turbine Number')
ylabel('Motion [m] or [deg]')
title(sprintf('FWF %s %s %s\n %s Fluctuations', farm_arrangement, farm_spacing, wave, DOF), 'interpreter', 'none')
ylim([-50, 50])



%% Tes with all turbines (fluctuations)


% full_tracking.(arrangement).(farm_spacing).(wave)(t).(DOF)waves

DOF = 'yaw_kal';
farm_arrangement = 'Staggered';
farm_spacing = 'SX50';
turbine = 9;
waves = fields(full_tracking.(farm_arrangement).(farm_spacing));

amplitudes = nan(1, length(waves));
for i = 1:length(waves)
    wave = waves{i};
    tmp = full_tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).(DOF);
    tmp = tmp - mean(tmp, 'all', 'omitnan');
    
    splts = split(wave, '_');

    if strcmp(splts{1}, 'LM0')
        amplitude = 0;
        wavelength = 0;
    else
        wavelength = wavelengths_from_string.(splts{1}) * 0.150;
        steepness = steepesses_from_string.(splts{2});
        amplitude = (steepness * wavelength) / (2 * pi);
    end
    amplitudes(i) = amplitude;
end



steepness_colors.AK06 = '#dd6e42';
steepness_colors.AK09 = '#e8dab2';
steepness_colors.AK12 = '#4f6d7a';


wavelength_markers.LM0 = 'o';
wavelength_markers.LM2 = '>';
wavelength_markers.LM3 = '<';
wavelength_markers.LM4 = 'p';
wavelength_markers.LM5 = 's';
wavelength_markers.LM25 = 'v';
wavelength_markers.LM33 = '^';

tickFontSize = 8;
labelFontSize = 10;



fit_x = nan(1, length(waves));
fit_y = nan(1, length(waves));

% Plot in order of increasing amplitude
[~, idxs] = sort(amplitudes, 'ascend');
reordered_waves = waves(idxs);
clc; close all
figure('color', 'white')
hold on
for i = 1:length(waves)
    wave = reordered_waves{i};
    tmp = full_tracking.(farm_arrangement).(farm_spacing).(wave)(turbine).(DOF);
    tmp = tmp - mean(tmp, 'all', 'omitnan');

    % crop to a certain number of STDs
    tmp(abs(tmp) > 3 * std(tmp, 0, 'all', 'omitnan')) = nan;
    
    splts = split(wave, '_');

    if strcmp(splts{1}, 'LM0')
        amplitude = 0;
        wavelength = 0;
        color = 'black';
    else
        wavelength = wavelengths_from_string.(splts{1}) * 0.150;
        steepness = steepesses_from_string.(splts{2});
        amplitude = (steepness * wavelength) / (2 * pi);

        color = steepness_colors.(splts{2});
    end

    VP = violinplot(amplitude * 100, abs(tmp), 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', 0.05, 'linewidth', 1);
    S = scatter(amplitude * 100, mean(abs(tmp), 'all', 'omitnan'), 50, wavelength_markers.(splts{1}), 'markerfacecolor', color, 'markeredgecolor', 'black', 'linewidth', 1);
    uistack(S, 'top')
    uistack(VP, 'bottom')

    fit_x(i) = amplitude * 100;
    fit_y(i) = mean(abs(tmp), 'all', 'omitnan');

end


P = polyfit(fit_x, fit_y, 2);
% plot_x = 0:0.1:max(fit_x);
plot_x = 0:0.1:2;
plot(plot_x, polyval(P, plot_x), 'color', 'black', 'linewidth', 2, 'linestyle', ':');
hold off


yline(0)
% xticks(1:12)
% xlabel('Turbine Number')
% ylabel('Motion [m] or [deg]')
title(sprintf('FWF %s %s %s\n %s Fluctuations', farm_arrangement, farm_spacing, wave, DOF), 'interpreter', 'none')
axis square
% yscale('log')
% ylim([-15, 15])
ylim([-2.5, 18])

ax = gca;
set(ax, 'FontSize', tickFontSize)
set(ax, 'TickLabelInterpreter', 'latex');
xlabel('$a$ [cm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(sprintf('$%s$ [%s]', DOF_symbs.(DOF), 'Deg'), 'interpreter', 'latex', 'fontsize', labelFontSize)
