%% Looking at intra-dynamics coupling for turbines
% Zein Sadek


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

% Experiment lengths
rotor_diameter = 0.15;


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
% TRY PLOTTING ONE DOF AGAINST ANOTHER (single case)
% SURGE + PITCH
% One spacing
% One wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spacing = 'SX50';
wave = 'LM5_AK06';

split_wave = split(wave, '_');
wavelength_name = split_wave{1};
steepness_name = split_wave{2};

steepness = steepesses_from_string.(steepness_name);
wavelength = wavelengths_from_string.(wavelength_name) * rotor_diameter;
amplitude = (wavelength * steepness) / (2 * pi);


clc; close all
figure('color', 'white')
tile = tiledlayout(1, length(centers));
sgtitle(sprintf('%s Floating Wind Farm:\n$\\lambda = %1.0fD, ak = %1.2f$', farm_arrangement, wavelength / rotor_diameter, steepness), 'Interpreter', 'latex')

for t = 1:length(centers)

    % Get signals
    turbine = centers(t);
    surge = tracking.(spacing).(wave)(turbine).x_kal;
    pitch = tracking.(spacing).(wave)(turbine).pitch_kal;
    
    % Remove mean
    surge = surge - mean(surge, 'all', 'omitnan');
    pitch = pitch - mean(pitch, 'all', 'omitnan');

    % Normalize by wave input
    normalized_surge = surge / amplitude;
    normalized_pitch = deg2rad(pitch) / steepness;
    
    % Plot
    h(t) = nexttile;
    colors = row_colors.(sprintf('Row%1.0f', t));

    skip = 1;
    hold on
    scatter(normalized_surge(1:skip:end), normalized_pitch(1:skip:end), 20, 'o', 'filled', ...
            'MarkerFaceColor', colors(1,:), 'MarkerFaceAlpha', 0.5)
    title(sprintf('Row %1.0f', t), 'interpreter', 'latex')
    hold off
    axis square


end

linkaxes(h, 'xy')
xlim([-3, 3])
ylim([-4, 4])

xlabel(tile, 'Normalized Surge: $x \mathbin{/} a$', 'interpreter', 'latex')
ylabel(tile, 'Normalized Pitch: $\theta \mathbin{/} ak$', 'interpreter', 'latex')






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRY PLOTTING ONE DOF AGAINST ANOTHER (multiple cases)
% SURGE + PITCH
% All spacings
% One wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wave = 'LM5_AK12';

split_wave = split(wave, '_');
wavelength_name = split_wave{1};
steepness_name = split_wave{2};

steepness = steepesses_from_string.(steepness_name);
wavelength = wavelengths_from_string.(wavelength_name) * rotor_diameter;
amplitude = (wavelength * steepness) / (2 * pi);

skip = 1;
marker_alpha = 0.1;

clc; close all
figure('color', 'white')
tile = tiledlayout(1, length(centers));
sgtitle(sprintf('%s Floating Wind Farm:\n$\\lambda = %1.0fD, ak = %1.2f$', farm_arrangement, wavelength / rotor_diameter, steepness), 'Interpreter', 'latex')

% Loop through rows of farm
for t = 1:length(centers)

    turbine = centers(t);
    colors  = row_colors.(sprintf('Row%1.0f', t));

    % Collect data over all spacings (pooled)
    normalized_surge_all = [];
    normalized_pitch_all = [];

    h(t) = nexttile;
    hold on

    % Loop through farm spacings (plot + pool)
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

        surge = tracking.(farm_spacing).(wave)(turbine).x_kal;
        pitch = tracking.(farm_spacing).(wave)(turbine).pitch_kal;

        % Remove mean
        surge = surge - mean(surge, 'omitnan');
        pitch = pitch - mean(pitch, 'omitnan');

        % Normalize by wave input
        normalized_surge = surge / amplitude;
        normalized_pitch = deg2rad(pitch) / steepness;

        % Pool (IMPORTANT: keep same decimation as plotted, if you want)
        xs = normalized_surge(1:skip:end);
        ys = normalized_pitch(1:skip:end);
        good = isfinite(xs) & isfinite(ys);

        normalized_surge_all = [normalized_surge_all; xs(good)];
        normalized_pitch_all = [normalized_pitch_all; ys(good)];

        % Plot
        scatter(xs, ys, 20, 'o', 'filled', ...
            'MarkerFaceColor', colors(s,:), 'MarkerFaceAlpha', marker_alpha, ...
            'HandleVisibility', 'off')

        % Legend proxy
        scatter(nan, nan, 20, 'o', 'filled', ...
            'MarkerFaceColor', colors(s,:), ...
            'DisplayName', sprintf('$S_x = %sD$', num2str(farm_spacings(s))))
    end

    title(sprintf('Row %1.0f', t), 'Interpreter', 'latex')
    axis square
    legend('Interpreter', 'latex', 'box', 'off', 'FontSize', 6)
    hold off
end


linkaxes(h, 'xy')
xlim([-3, 3])
ylim([-4, 4])

xlabel(tile, 'Normalized Surge: $x \mathbin{/} a$', 'interpreter', 'latex')
ylabel(tile, 'Normalized Pitch: $\theta \mathbin{/} ak$', 'interpreter', 'latex')







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRY PLOTTING ONE DOF AGAINST ANOTHER (multiple cases)
% PCA Overlayed
% SURGE + PITCH
% All spacings
% One wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wave = 'LM5_AK12';

split_wave = split(wave, '_');
wavelength_name = split_wave{1};
steepness_name = split_wave{2};

steepness = steepesses_from_string.(steepness_name);
wavelength = wavelengths_from_string.(wavelength_name) * rotor_diameter;
amplitude = (wavelength * steepness) / (2 * pi);

skip = 1;
marker_alpha = 0.1;

clc; close all
figure('color', 'white')
tile = tiledlayout(1, length(centers));
sgtitle(sprintf('%s Floating Wind Farm:\n$\\lambda = %1.0fD, ak = %1.2f$', farm_arrangement, wavelength / rotor_diameter, steepness), 'Interpreter', 'latex')

% Loop through rows of farm
for t = 1:length(centers)

    turbine = centers(t);
    colors  = row_colors.(sprintf('Row%1.0f', t));

    % Collect data over all spacings (pooled)
    normalized_surge_all = [];
    normalized_pitch_all = [];

    h(t) = nexttile;
    hold on

    % Loop through farm spacings (plot + pool)
    for s = 1:length(farm_spacings)
        farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

        surge = tracking.(farm_spacing).(wave)(turbine).x_kal;
        pitch = tracking.(farm_spacing).(wave)(turbine).pitch_kal;

        % Remove mean
        surge = surge - mean(surge, 'omitnan');
        pitch = pitch - mean(pitch, 'omitnan');

        % Normalize by wave input
        normalized_surge = surge / amplitude;
        normalized_pitch = deg2rad(pitch) / steepness;

        % Pool (IMPORTANT: keep same decimation as plotted, if you want)
        xs = normalized_surge(1:skip:end);
        ys = normalized_pitch(1:skip:end);
        good = isfinite(xs) & isfinite(ys);

        % Edge case fix
        if length(good) < 3600
            good = ones(1,3600);
            xs(3600) = xs(end);
            ys(3600) = ys(end);
            disp('Messed up')
        end

        normalized_surge_all = [normalized_surge_all; xs(good)];
        normalized_pitch_all = [normalized_pitch_all; ys(good)];

        % Plot
        scatter(xs, ys, 20, 'o', 'filled', ...
            'MarkerFaceColor', colors(s,:), 'MarkerFaceAlpha', marker_alpha, ...
            'HandleVisibility', 'off')

        % Legend proxy
        scatter(nan, nan, 20, 'o', 'filled', ...
            'MarkerFaceColor', colors(s,:), ...
            'DisplayName', sprintf('$S_x = %sD$', num2str(farm_spacings(s))))
    end

    % ---- PCA on pooled cloud (per row) ----
    A = pca_principal_axes_2d(normalized_surge_all(:), normalized_pitch_all(:));

    % Overlay principal axes + ellipse
    plot_pca_axes_ellipse(A, 'k', 1);

    title(sprintf('Row %1.0f', t), 'Interpreter', 'latex')
    axis square
    legend('Interpreter', 'latex', 'box', 'off', 'FontSize', 6)
    hold off
end


linkaxes(h, 'xy')
xlim([-3, 3])
ylim([-4, 4])

xlabel(tile, 'Normalized Surge: $x \mathbin{/} a$', 'interpreter', 'latex')
ylabel(tile, 'Normalized Pitch: $\theta \mathbin{/} ak$', 'interpreter', 'latex')







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE PCA FOR EACH WAVE, OVER ALL SPACINGS
% SURGE + PITCH
% All spacings
% One wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
% Loop through wave steepnesses
for st = 1:length(wave_steepnesses)

    % Loop through wavelengths
    for w = 1:length(wavelengths)
        steepness = wave_steepnesses(st);
        steep = compose('%02d', round(100 * steepness));
        wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

        wavelength = wavelengths(w) * rotor_diameter;
        amplitude = (wavelength * steepness) / (2 * pi);
        
        skip = 1;
               
        % Loop through rows of farm
        for t = 1:length(centers)
        
            turbine = centers(t);
            colors  = row_colors.(sprintf('Row%1.0f', t));
        
            % Collect data over all spacings (pooled)
            normalized_surge_all = [];
            normalized_pitch_all = [];
        
            % Loop through farm spacings (plot + pool)
            for s = 1:length(farm_spacings)
                farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
        
                surge = tracking.(farm_spacing).(wave)(turbine).x_kal;
                pitch = tracking.(farm_spacing).(wave)(turbine).pitch_kal;
        
                % Remove mean
                surge = surge - mean(surge, 'omitnan');
                pitch = pitch - mean(pitch, 'omitnan');
        
                % Normalize by wave input
                normalized_surge = surge / amplitude;
                normalized_pitch = deg2rad(pitch) / steepness;
        
                % Pool (IMPORTANT: keep same decimation as plotted, if you want)
                xs = normalized_surge(1:skip:end);
                ys = normalized_pitch(1:skip:end);
                good = isfinite(xs) & isfinite(ys);
        
                % Edge case fix
                if length(good) < 3600
                    good = ones(1,3600);
                    xs(3600) = xs(end);
                    ys(3600) = ys(end);
                    disp('Messed up')
                end
        
                normalized_surge_all = [normalized_surge_all; xs(good)];
                normalized_pitch_all = [normalized_pitch_all; ys(good)];
            end
        
            % ---- PCA on pooled cloud (per row) ----
            A = pca_principal_axes_2d(normalized_surge_all(:), normalized_pitch_all(:));

            % Save PCA results
            PrincipalAxes.(wave)(t) = A;
        
        end

    end
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AXIS RATIO VS WAVELENGTH
% SURGE + PITCH
% All spacings
% One wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('color', 'white')
tile = tiledlayout(1,4);
sgtitle(sprintf('%s Floating Wind Farm: Principal Axes Ratio\nSurge-Pitch', farm_arrangement), 'Interpreter', 'latex')
% Loop through rows of turbines
for t = 1:4
    turbine = centers(t);

    h(t) = nexttile;
    hold on

    % Loop through wavelengths
    tmptmp = nan(1,4);
    for w = 1:length(wavelengths)
        % Loop through wave steepnesses
        tmp = nan(1,3);
        for st = 1:length(wave_steepnesses)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            if st == 1
                marker = '^';
            elseif st == 2
                marker = 'square';
            elseif st == 3
                marker = 'o';
            end

            % Plot
            scatter(wavelengths(w), PrincipalAxes.(wave)(t).axis_ratio, ...
                    marker, 'filled', 'HandleVisibility', 'off')

            tmp(st) = PrincipalAxes.(wave)(t).axis_ratio;
        end
        scatter(wavelengths(w), mean(tmp, 'all', 'omitnan'), ...
                    marker, 'filled', 'markerfacecolor', 'k', ...
                    'HandleVisibility', 'off')
        tmptmp(w) = mean(tmp, 'all', 'omitnan');
    end
    plot(wavelengths, tmptmp, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off')
    hold off
    title(sprintf('Row %1.0f', t), 'Interpreter', 'latex')
    axis square
end

% Legend
hold on
for st = 1:3
    if st == 1
        marker = '^';
    elseif st == 2
        marker = 'square';
    elseif st == 3
        marker = 'o';
    end

    steepness = wave_steepnesses(st);
    label = sprintf('$ak = %1.2f$', steepness);
    scatter(nan, nan, marker, 'filled', 'MarkerFaceColor', 'k', 'displayname', label)
end

leg = legend('box', 'off', 'interpreter', 'latex');
leg.Layout.Tile = 'east';
linkaxes(h, 'xy')
ylim([0,1])
xlim([1.5, 5.5])

xlabel(tile, '$\lambda \mathbin{/} D$', 'interpreter', 'latex')
ylabel(tile, '$\sigma_{2} \mathbin{/} \sigma_{1}$', 'interpreter', 'latex')







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING AXIS ANGLE VS WAVELENGTH
% SURGE + PITCH
% All spacings
% One wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('color', 'white')
tile = tiledlayout(1,4);
sgtitle(sprintf('%s Floating Wind Farm: Principal Axes Angle\nSurge-Pitch', farm_arrangement), 'Interpreter', 'latex')
% Loop through rows of turbines
for t = 1:4
    turbine = centers(t);
    h(t) = nexttile;
    hold on
    

    % Loop through wavelengths
    tmptmp = nan(1,4);
    for w = 1:length(wavelengths)
        tmp = nan(1,3);
        % Loop through wave steepnesses
        for st = 1:length(wave_steepnesses)
        
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            if st == 1
                marker = '^';
            elseif st == 2
                marker = 'square';
            elseif st == 3
                marker = 'o';
            end

            % Plot
            angle = PrincipalAxes.(wave)(t).angle_deg;
            if angle < 0
                angle = angle + 360;
            end
            scatter(wavelengths(w), angle, ...
                    marker, 'filled', 'HandleVisibility', 'off')
            tmp(st) = angle;
        end
        scatter(wavelengths(w), mean(tmp, 'all', 'omitnan'), ...
                    marker, 'filled', 'markerfacecolor', 'k', ...
                    'HandleVisibility', 'off')
        tmptmp(w) = mean(tmp, 'all', 'omitnan');
    end
    plot(wavelengths, tmptmp, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off')
    hold off
    title(sprintf('Row %1.0f', t), 'Interpreter', 'latex')
    axis square
end

% Legend
hold on
for st = 1:3
    if st == 1
        marker = '^';
    elseif st == 2
        marker = 'square';
    elseif st == 3
        marker = 'o';
    end

    steepness = wave_steepnesses(st);
    label = sprintf('$ak = %1.2f$', steepness);
    scatter(nan, nan, marker, 'filled', 'MarkerFaceColor', 'k', 'displayname', label)
end

leg = legend('box', 'off', 'interpreter', 'latex');
leg.Layout.Tile = 'east';
linkaxes(h, 'xy')
ylim([0,360])
xlim([1.5, 5.5])

xlabel(tile, '$\lambda \mathbin{/} D$', 'interpreter', 'latex')
ylabel(tile, '$\alpha$ [deg]', 'interpreter', 'latex')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING EXPLAINED VARIANCE VS WAVELENGTH
% SURGE + PITCH
% All spacings
% One wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('color', 'white')
tile = tiledlayout(1,4);
sgtitle(sprintf('%s Floating Wind Farm: Principal Axes Explained Variance', farm_arrangement), 'Interpreter', 'latex')
% Loop through rows of turbines
for t = 1:4
    turbine = centers(t);

    h(t) = nexttile;
    hold on

    % Loop through wavelengths
    tmptmp = nan(1,4);
    for w = 1:length(wavelengths)
        tmp = nan(1,3);
        % Loop through wave steepnesses
        for st = 1:length(wave_steepnesses)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            if st == 1
                marker = '^';
            elseif st == 2
                marker = 'square';
            elseif st == 3
                marker = 'o';
            end

            % Plot
            scatter(wavelengths(w), PrincipalAxes.(wave)(t).explained, ...
                    marker, 'filled', 'HandleVisibility', 'off')
            tmp(st) = PrincipalAxes.(wave)(t).explained;
        end
        scatter(wavelengths(w), mean(tmp, 'all', 'omitnan'), ...
                    marker, 'filled', 'markerfacecolor', 'k', ...
                    'HandleVisibility', 'off')
        tmptmp(w) = mean(tmp, 'all', 'omitnan');
    end
    plot(wavelengths, tmptmp, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off')
    title(sprintf('Row %1.0f', t), 'Interpreter', 'latex')
end

% Legend
hold on
for st = 1:3
    if st == 1
        marker = '^';
    elseif st == 2
        marker = 'square';
    elseif st == 3
        marker = 'o';
    end

    steepness = wave_steepnesses(st);
    label = sprintf('$ak = %1.2f$', steepness);
    scatter(nan, nan, marker, 'filled', 'MarkerFaceColor', 'k', 'displayname', label)
end

leg = legend('box', 'off', 'interpreter', 'latex');
leg.Layout.Tile = 'east';
linkaxes(h, 'xy')
ylim([0,1])
xlim([1.5, 5.5])

xlabel(tile, '$\lambda \mathbin{/} D$', 'interpreter', 'latex')
% ylabel(tile, '$\alpha$ [deg]', 'interpreter', 'latex')










%% Functions



function A = pca_principal_axes_2d(x, y)

    conf = 0.90;        % confidence for ellipse
    use_core = true;    % set false if you want full PCA
    core_frac = 0.90;   % central fraction to keep

    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);

    X = [x y];

    % ----- Core selection (robust against triangle clouds) -----
    if use_core
        mu0 = median(X,1);
        C0  = cov(X);
        d2  = mahal(X, X);
        thr = quantile(d2, core_frac);
        X   = X(d2 <= thr,:);
    end

    mu = mean(X,1)';
    C  = cov(X);

    [V,D] = eig(C);
    [lam, idx] = sort(diag(D), 'descend');
    V = V(:,idx);

    % Ellipse scaling (2 DOF)
    s = sqrt(chi2inv(conf,2));

    A.mu = mu;
    A.V  = V;
    A.lam = lam;
    A.axis_sigma = sqrt(lam);
    A.axis_len   = s * sqrt(lam);
    A.angle_deg  = atan2d(V(2,1), V(1,1));
    A.axis_ratio = A.axis_sigma(2)/A.axis_sigma(1);
    A.explained  = lam(1)/sum(lam);
end



function d2 = mahal2(X, mu, C)
% squared Mahalanobis distance for each row of X (2D)
    Xc = X - mu;
    C = (C + C')/2;
    [R,p] = chol(C);
    if p ~= 0
        C = C + 1e-8*eye(2);
        R = chol(C);
    end
    Z = Xc / R;
    d2 = sum(Z.^2,2);
end

function plot_pca_axes_ellipse(A, lineColor, lineWidth)
% Plots PCA major/minor axes and confidence ellipse.

    if nargin < 2, lineColor = 'k'; end
    if nargin < 3, lineWidth = 2; end

    mu = A.mu;
    V  = A.V;
    L  = A.axis_len;

    % Axes lines
    plot([mu(1)-L(1)*V(1,1), mu(1)+L(1)*V(1,1)], ...
         [mu(2)-L(1)*V(2,1), mu(2)+L(1)*V(2,1)], ...
         ':', 'Color', lineColor, 'LineWidth', lineWidth, 'HandleVisibility','off');

    plot([mu(1)-L(2)*V(1,2), mu(1)+L(2)*V(1,2)], ...
         [mu(2)-L(2)*V(2,2), mu(2)+L(2)*V(2,2)], ...
         ':', 'Color', lineColor, 'LineWidth', lineWidth, 'HandleVisibility','off');

    % Ellipse
    t = linspace(0, 2*pi, 400);
    E = (V * diag(L) * [cos(t); sin(t)]) + mu;
    plot(E(1,:), E(2,:), '-', 'Color', lineColor, 'LineWidth', lineWidth-0.5, ...
         'HandleVisibility','off');
end

