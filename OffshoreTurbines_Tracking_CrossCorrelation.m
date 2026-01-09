%% Looking at cross-correlation of motion
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
% TEST CROSS CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test cross-correlation
farm_spacing = 'SX50';
wave = 'LM5_AK12';
DOF = 'pitch_kal';

row1_signal = tracking.(farm_spacing).(wave)(2).(DOF);
row2_signal = tracking.(farm_spacing).(wave)(5).(DOF);


clc;
fs = 30;          % Hz
maxLag_s = 120;    % seconds
out = xcorr_metrics(row1_signal, row2_signal, fs, maxLag_s);
disp(out)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG MAX CORRELATION VALUE
% FOR ALL SPACINGS WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc;
% fs = 30;          % Hz
% maxLag_s = 10;    % seconds
% 
% % Specifying the 'fixed' turbine row and show the max correlation value
% % with all other waked rows
% centers = [2, 5, 8, 12];
% fixed_row = 1;
% DOF = 'pitch_kal';
% 
% % Determine waked turbines
% waked_turbines = centers(fixed_row + 1:end);
% num_waked_turbines = length(waked_turbines);
% 
% % Plot things
% wave_steepnesses = [0.06, 0.09, 0.12];
% wavelengths = [5,4,3,2];
% steepness_alpha = [0.3, 0.6, 1];
% sz = 100;
% spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
% wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};
% 
% 
% clc; close all
% figure('color', 'white')
% t = tiledlayout(1, num_waked_turbines);
% sgtitle(sprintf('Correlation coefficient of %s ($%s$)', names.(DOF), symbs.(DOF)), 'interpreter', 'latex')
% 
% % Loop through waked turbines
% for c = 1:num_waked_turbines
% 
%     reference_turbine = centers(fixed_row);
%     correlating_turbine = waked_turbines(c);
% 
%     h(c) = nexttile;
%     title(sprintf('Row %1.0f to Row %1.0f', fixed_row, ceil(correlating_turbine / 3)))
%     hold on
%     % Loop through spacings
%     for s = 1:length(farm_spacings)
%         farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
%         caze = strcat("WT60_", farm_spacing, "_AG0");
%         fprintf('%s\n', caze)
% 
%         tmp = tracking.(farm_spacing);
%         waves = fieldnames(tmp);
% 
%         % Loop through waves
%         for st = 1:length(wave_steepnesses)
%             wave_steepness = wave_steepnesses(st);
%             steep = compose('%02d', round(100 * wave_steepness));
%             disp(steep{1})
% 
%             for w = 1:length(wavelengths)
%                 % Make wave
%                 wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
%                 harmonic_ratio = farm_spacings(s) / wavelengths(w);
% 
% 
%                 % Get signals
%                 turbine_A_signal = tracking.(farm_spacing).(wave)(reference_turbine).(DOF);
%                 turbine_B_signal = tracking.(farm_spacing).(wave)(correlating_turbine).(DOF);
% 
% 
%                 % Make the same length
%                 n = min(numel(turbine_A_signal), numel(turbine_B_signal));
%                 A = turbine_A_signal(1:n);
%                 B = turbine_B_signal(1:n);
% 
% 
%                 % Cross-correlate
%                 out = xcorr_metrics(A, B, fs, maxLag_s);
% 
%                 % out.rho_max ~ largest magnitude XC coefficient
%                 % out.rho_abs_max ~ abs of largest magnitude XC coefficient
%                 % out.tau_max ~ time lag at largest peak
% 
%                 scatter(harmonic_ratio, out.rho_max, sz, spacing_shapes{s}, 'filled', ...
%                         'MarkerFaceColor', wave_colors{w}, 'markerfacealpha', steepness_alpha(st))
% 
%             end
%         end
%     end
%     xline(0.5:1:2.5, 'linestyle', '--')
%     xline(1:1:2, 'linestyle', '-')
%     hold off
% end
% 
% linkaxes(h, 'xy')
% xlim([0.4, 2.5])
% xlabel(t, '$S_x / \lambda$', 'interpreter', 'latex')
% ylabel(t, sprintf('$\\rho_{%s}$', symbs.(DOF)), 'interpreter', 'latex')






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTIG ABSOLUTE VALUE OF MAX CORRELATION VALUE
% FOR ALL SPACINGS WAVES
% SPECIFYING FIXED ROW AND DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
fs = 30;          % Hz
maxLag_s = 100;   % seconds

% Specifying the 'fixed' turbine row and show the max correlation value
% with all other waked rows
% centers = [2, 5, 8, 12];
if strcmp(farm_arrangement, 'Inline')
    centers = [2,5,8,12];
else
    centers = [2,4,7,9];
end
fixed_row = 2;
DOF = 'pitch_kal';

% Determine waked turbines
waked_turbines = centers(fixed_row + 1:end);
num_waked_turbines = length(waked_turbines);

% Plot things
wave_steepnesses = [0.06, 0.09, 0.12];
wavelengths = [5,4,3,2];
steepness_alpha = [0.3, 0.6, 1];
sz = 100;
spacing_shapes = {'o', 'diamond', '^', 'v', 'square'};
wave_colors = {'#EC4E20', '#FF9505', '#4C4B63', '#ABA8B2'};


clc; close all
figure('color', 'white')
t = tiledlayout(1, num_waked_turbines);
sgtitle(sprintf('ABS Correlation coefficient of %s ($%s$): $|\\rho_{%s}|$', names.(DOF), symbs.(DOF), symbs.(DOF)), 'interpreter', 'latex')

% Loop through waked turbines
for c = 1:num_waked_turbines

    reference_turbine = centers(fixed_row);
    correlating_turbine = waked_turbines(c);

    h(c) = nexttile;
    title(sprintf('Row %1.0f to Row %1.0f', fixed_row, ceil(correlating_turbine / 3)))
    hold on
    % Loop through spacings
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        caze = strcat("WT60_", farm_spacing, "_AG0");
        fprintf('%s\n', caze)
        
        tmp = tracking.(farm_spacing);
        waves = fieldnames(tmp);

        % Loop through waves
        for st = 1:length(wave_steepnesses)
            wave_steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * wave_steepness));
            disp(steep{1})

            for w = 1:length(wavelengths)
                % Make wave
                wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
                harmonic_ratio = farm_spacings(s) / wavelengths(w);


                % Get signals
                turbine_A_signal = tracking.(farm_spacing).(wave)(reference_turbine).(DOF);
                turbine_B_signal = tracking.(farm_spacing).(wave)(correlating_turbine).(DOF);
                
    
                % Make the same length
                n = min(numel(turbine_A_signal), numel(turbine_B_signal));
                A = turbine_A_signal(1:n);
                B = turbine_B_signal(1:n);
    
    
                % Cross-correlate
                out = xcorr_metrics(A, B, fs, maxLag_s);
    
                % out.rho_max ~ largest magnitude XC coefficient
                % out.rho_abs_max ~ abs of largest magnitude XC coefficient
                % out.tau_max ~ time lag at largest peak

                scatter(harmonic_ratio, out.rho_abs_max, sz, spacing_shapes{s}, 'filled', ...
                        'MarkerFaceColor', wave_colors{w}, 'markerfacealpha', steepness_alpha(st), ...
                        'HandleVisibility', 'off')

            end
        end
    end

     % %%% Legend
    if c == 1
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

linkaxes(h, 'xy')
xlim([0.5, 2.5])
ylim([0, 1])
xlabel(t, '$S_x / \lambda$', 'interpreter', 'latex')
ylabel(t, sprintf('$| \\rho_{%s} |$', symbs.(DOF)), 'interpreter', 'latex')

















%% Functions

function out = xcorr_metrics(x, y, fs, maxLag_s)
% Returns peak correlation and lag at peak (seconds)

x = x(:); y = y(:);
x = x - mean(x, 'omitnan');
y = y - mean(y, 'omitnan');

% Optional: handle NaNs by simple fill (or remove segments)
x = fillmissing(x,'linear','EndValues','nearest');
y = fillmissing(y,'linear','EndValues','nearest');

maxLag = round(maxLag_s * fs);
[r, lags] = xcorr(x, y, maxLag, 'coeff');
tau = lags / fs;

% Peak magnitude (use abs so anti-correlation counts as strong coupling)
[~, idx] = max(abs(r));
out.rho_max = r(idx);
out.rho_abs_max = abs(r(idx));
out.tau_max = tau(idx);

% Also useful: zero-lag correlation
[~, iz] = min(abs(tau));
out.rho0 = r(iz);
end
