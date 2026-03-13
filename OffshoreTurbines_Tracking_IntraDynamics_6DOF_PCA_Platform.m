%% Looking at intra-dynamics coupling for turbines
% Zein Sadek


clear; close all; clc;
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM")
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")
addpath('/Users/zeinsadek/Documents/MATLAB/SpiderPlot')

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

% Platform length (x)
floater_spacing = 0.1;
platform_length_x = floater_spacing * cosd(30);
hub_height_from_platform = 0.1;
hub_height_from_surface = 0.117;


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
% 6 DOF PCA TO SEE WHICH DOFS ARE CORRELATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

turbine = centers(3);
farm_spacing = 'SX50';
wave = 'LM5_AK12';

time = tracking.(farm_spacing).(wave)(turbine).time;
% Translations
surge = tracking.(farm_spacing).(wave)(turbine).x_kal;
heave = tracking.(farm_spacing).(wave)(turbine).y_kal;
sway  = tracking.(farm_spacing).(wave)(turbine).z_kal;

% Rotations
roll  = tracking.(farm_spacing).(wave)(turbine).roll_kal;
pitch = tracking.(farm_spacing).(wave)(turbine).pitch_kal;
yaw   = tracking.(farm_spacing).(wave)(turbine).yaw_kal;


% Normalize
% Translations by platform length
surge_norm = surge / floater_spacing;
heave_norm = heave / floater_spacing;
sway_norm  = sway  / floater_spacing;


% Normalize rotations by displacement at hub, over platform length
roll_norm  = (deg2rad(roll) * hub_height_from_platform) / floater_spacing;
pitch_norm = (deg2rad(pitch) * hub_height_from_platform) / floater_spacing;
yaw_norm   = deg2rad(yaw) / sqrt(3);


% Plot of signals
figure('color', 'white')
tile = tiledlayout(6,1);
sgtitle('Platform Scaling')

h(1) = nexttile;
plot(time, surge_norm)
ax = gca;
ax.XAxis.Visible = 'off';
box off
title('Surge')

h(2) = nexttile;
plot(time, heave_norm - mean(heave_norm))
ax = gca;
ax.XAxis.Visible = 'off';
box off
title('Heave')

h(3) = nexttile;
plot(time, sway_norm - mean(sway_norm))
ax = gca;
ax.XAxis.Visible = 'off';
box off
title('Sway')

h(4) = nexttile;
plot(time, roll_norm - mean(roll_norm))
ax = gca;
ax.XAxis.Visible = 'off';
box off
title('Roll')

h(5) = nexttile;
plot(time, pitch_norm - mean(pitch_norm))
ax = gca;
ax.XAxis.Visible = 'off';
box off
title('Pitch')

h(6) = nexttile;
plot(time, yaw_norm - mean(yaw_norm))
box off
title('Yaw')

linkaxes(h, 'xy')
ylim([-0.25, 0.25])
xlabel(tile, 'Time [s]')
ylabel(tile, 'Displacement [m]')


%%

% Assemble data matrix
X = [surge_norm(:), heave_norm(:), sway_norm(:), ...
     roll_norm(:),  pitch_norm(:), yaw_norm(:)];
X = X - mean(X, 1, 'omitnan');


% ------------------------------------------------------------
% PLOTS: explained variance, loadings, correlation, modal coords
% ------------------------------------------------------------

dof_names = {'surge','heave','sway','roll','pitch','yaw'};

% 0) (Recommended) Remove rows with NaNs before plotting / cov
ok = all(isfinite(X),2);
Xc = X(ok,:);                     % cleaned X
Cc = cov(Xc);                     % cleaned covariance
Rc = corrcoef(Xc);                % correlation matrix

% Recompute eig on cleaned covariance (optional but recommended)
[V,D] = eig(Cc);
lambda = diag(D);
[lambda,idx] = sort(lambda,'descend');
V = V(:,idx);
explained = lambda/sum(lambda);

% 1) Explained variance bar plot
figure('Color','w'); 
bar(100*explained,'FaceAlpha',0.9);
grid on; box on;
xlabel('Principal component');
ylabel('Explained variance (%)');
title(sprintf('Platform Scaling\n6-DOF PCA explained variance (%s, %s, turbine %d)', farm_spacing, wave, turbine), 'interpreter', 'none');
xlim([0.5 6.5]);
ylim([0, 100])

figure('color', 'white')
plot(100*cumsum(explained), 1:6, '-o', 'LineWidth', 1.5);



% 2) Mode loading plots (PC1–PC3)
nModesToShow = 3;

figure('Color','w');
tiledlayout(3,1)
for k = 1:nModesToShow
    h(k) = nexttile;
    stem(abs(V(:,k)),'filled','LineWidth',1.2);
    grid on; box on;
    ylabel(sprintf('PC%d',k));
    ylim([-1 1]);
    set(gca,'XTick',1:6,'XTickLabel',dof_names);
    if k==1
        title(sprintf('PCA mode loadings (eigenvectors) — %s, %s', farm_spacing, wave));
    end
    % annotate explained variance
    text(0.98,0.85,sprintf('%.1f%%',100*explained(k)), ...
        'Units','normalized','HorizontalAlignment','right','FontWeight','bold');
    ylim([0, 1.2])
end
xlabel('DOF');
linkaxes(h, 'xy')


% 4) Heatmap: correlation matrix
figure('Color','w');
imagesc(corr(X));
axis equal tight;
colorbar;
set(gca,'XTick',1:6,'XTickLabel',dof_names,'YTick',1:6,'YTickLabel',dof_names);
title('Correlation matrix corr(X)');
% Optional numbers on cells:
for i=1:6
    for j=1:6
        text(j,i,sprintf('%.2f',Rc(i,j)),'HorizontalAlignment','center');
    end
end

% 6) (Optional, very useful) Contribution "weights" per mode
% Using squared loadings as percent contribution inside each mode:
W = V(:,1:nModesToShow).^2;
W = W ./ sum(W,1);  % each column sums to 1

figure('Color','w');
bar(100*W','stacked');
grid on; box on;
set(gca,'XTick',1:nModesToShow,'XTickLabel',compose('PC%d',1:nModesToShow));
ylabel('Within-mode contribution (%)');
legend(dof_names,'Location','eastoutside');
title(sprintf('Platform Scaling\nWithin-mode DOF contributions (squared loadings)'));







%% Loop through cases and plot useful metrics


turbine = centers(3);

% Loop through farm spacings
for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    % Loop through wave steepnesses
    for st = 1:length(wave_steepnesses)
    
        % Loop through wavelengths
        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            wavelength = wavelengths(w) * rotor_diameter;
            amplitude = (wavelength * steepness) / (2 * pi);

            % Translations
            surge = tracking.(farm_spacing).(wave)(turbine).x_kal;
            heave = tracking.(farm_spacing).(wave)(turbine).y_kal;
            sway  = tracking.(farm_spacing).(wave)(turbine).z_kal;
            
            % Rotations
            roll  = tracking.(farm_spacing).(wave)(turbine).roll_kal;
            pitch = tracking.(farm_spacing).(wave)(turbine).pitch_kal;
            yaw   = tracking.(farm_spacing).(wave)(turbine).yaw_kal;
            
            
            % % Normalize: Platform
            % Translations by platform length
            surge_norm = surge / floater_spacing;
            heave_norm = heave / floater_spacing;
            sway_norm  = sway  / floater_spacing;

            % Normalize rotations by displacement at hub, over platform length
            roll_norm  = (deg2rad(roll) * hub_height_from_surface) / floater_spacing;
            pitch_norm = (deg2rad(pitch) * hub_height_from_surface) / floater_spacing;
            yaw_norm   = deg2rad(yaw) / sqrt(3);

            % Normalize: Wave
            % Translations by platform length
            % surge_norm = surge / amplitude;
            % heave_norm = heave / amplitude;
            % sway_norm  = sway  / amplitude;
            % 
            % % Normalize rotations by displacement at hub, over platform length
            % roll_norm  = deg2rad(roll) / steepness;
            % pitch_norm = deg2rad(pitch) / steepness;
            % yaw_norm   = deg2rad(yaw) / steepness;
            
            
            % Assemble data matrix
            X = [surge_norm(:), heave_norm(:), sway_norm(:), ...
                 roll_norm(:),  pitch_norm(:), yaw_norm(:)];
            X = X - mean(X, 1, 'omitnan');


            % 0) (Recommended) Remove rows with NaNs before plotting / cov
            ok = all(isfinite(X),2);
            Xc = X(ok,:);                     % cleaned X
            Cc = cov(Xc);                     % cleaned covariance
            Rc = corrcoef(Xc);                % correlation matrix
            
            % Recompute eig on cleaned covariance (optional but recommended)
            [V,D] = eig(Cc);
            lambda = diag(D);
            [lambda,idx] = sort(lambda,'descend');
            V = V(:,idx);
            explained = lambda/sum(lambda);

            % Save
            decomposition.(farm_spacing).(wave).eigenvalues = lambda;
            decomposition.(farm_spacing).(wave).eigenvectors = V;
            decomposition.(farm_spacing).(wave).fractional_energy = explained;

        end
    end
end


%% Plot energy per mode for all cases


wave_colors = slanCM('parula', 5);
steepness_linestyle = {':', '--', '-'};

clc; close all
figure('color', 'white')
tile = tiledlayout(1,1);

% Loop through farm spacings
% for s = 1:length(farm_spacings)
for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    h(s) = nexttile;
    hold on
    % Loop through wave steepnesses
    % for st = 1:length(wave_steepnesses)
    for st = 1:3
    
        % Loop through wavelengths
        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            % plot(100 * cumsum(decomposition.(farm_spacing).(wave).fractional_energy), 1:6, ...
            %      'linewidth', 2, 'color', wave_colors(w,:), 'linestyle', steepness_linestyle{st}, ...
            %      'HandleVisibility', 'off')

            plot(1:6, decomposition.(farm_spacing).(wave).fractional_energy, ...
                 'linewidth', 2, 'color', wave_colors(w,:), 'linestyle', steepness_linestyle{st}, ...
                 'HandleVisibility', 'off')

        end
    end
    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
    xticks(1:6)
end

hold on
% Legend for wave colors
for w = 1:length(wavelengths)
    label = sprintf('$\\lambda = %sD$', num2str(wavelengths(w)));
    plot(nan, nan, 'linewidth', 2, 'color', wave_colors(w,:), 'displayname', label)
end

% Legend for line style
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
for st = 1:3
    label = sprintf('$ak = %0.2f$', wave_steepnesses(st));
    plot(nan, nan, 'linewidth', 2, 'color', 'k', 'displayname', label, 'linestyle', steepness_linestyle{st})
end
hold off

leg = legend('box', 'off', 'interpreter', 'latex', 'orientation', 'vertical', 'location', 'northeast');
% leg.Layout.Tile = 'east';

% yscale('log')
ylim([0,1])
linkaxes(h, 'xy')
xlabel(tile, 'Modes', 'interpreter', 'latex')
ylabel(tile, 'Captured Variance [\%]', 'interpreter', 'latex')



%% Look at the energy contained within the first mode vs wave parameters


mode = 1;

clc; close all
figure('color', 'white')
tile = tiledlayout(1,1);
sgtitle(sprintf('Mode %1.0f', mode))

% Loop through farm spacings
% for s = 1:length(farm_spacings)
for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    h(s) = nexttile;
    hold on
    % Loop through wave steepnesses
    % for st = 1:length(wave_steepnesses)
    for st = 1:3
    
        tmpx = nan(1,length(wavelengths));
        tmpy = nan(1,length(wavelengths));
        % Loop through wavelengths
        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            wavelength = wavelengths(w) * rotor_diameter;
            wavenumber = (2 * pi) / (wavelength);
            amplitude = (wavelength * steepness) / (2 * pi);

            %%% Different geometric scalings
            % Phase difference between front/rear buoys
            % buoy_phase = rad2deg((2 * pi * floater_spacing) / (wavelength));
            buoy_phase = (floater_spacing * cosd(30)) / wavelength;

            % Max geometric pitch angle
            geometric_pitch = rad2deg((2 * amplitude / floater_spacing) * sin((wavenumber * floater_spacing) / 2));


            % Plot
            scatter(buoy_phase, decomposition.(farm_spacing).(wave).fractional_energy(mode), ...
                    50, 'filled', 'markerfacecolor', wave_colors(w,:),...
                    'HandleVisibility', 'off')

            tmpx(w) = buoy_phase;
            tmpy(w) = decomposition.(farm_spacing).(wave).fractional_energy(mode);
        end

        % Plot connecting line
        [sorted_tmpx, idx] = sort(tmpx, 'ascend');
        P = plot(sorted_tmpx, tmpy(idx), 'color', 'k', ...
                 'linestyle', steepness_linestyle{st}, ...
                 'linewidth', 1.5, ...
                 'HandleVisibility', 'off');
        uistack(P, 'bottom')
    end
    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
end

% Legends
hold on
for w = 1:length(wavelengths)
    plot(nan,nan,'LineWidth',2,'Color',wave_colors(w,:), ...
        'DisplayName',sprintf('$\\lambda=%dD$', wavelengths(w)));
end
plot(nan,nan,'Color','w','DisplayName',' ')
for st = 1:length(wave_steepnesses)
    plot(nan,nan,'k','LineWidth',2,'LineStyle',steepness_linestyle{st}, ...
        'DisplayName',sprintf('$ak=%0.2f$', wave_steepnesses(st)));
end
hold off
leg = legend('Interpreter','latex','Box','off','location', 'northeast'); 
% leg.Layout.Tile = 'east';


linkaxes(h, 'xy')
ylim([0,1])
xlim([0,0.5])
ylabel(tile, 'Mode 1 Captured Varience [\%]', 'interpreter', 'latex')
xlabel(tile, 'Front/Rear Buoy Phase Misalingment: $s_x \mathbin{/} \lambda$', 'interpreter', 'latex')




%% Look at the maximum contribution from a single DOF in mode 1


clc; close all
figure('color', 'white')
tile = tiledlayout(1,1);

% Loop through farm spacings
% for s = 1:length(farm_spacings)
for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    h(s) = nexttile;
    hold on
    % Loop through wave steepnesses
    % for st = 1:length(wave_steepnesses)
    for st = 1:3
    
        tmpx = nan(1,length(wavelengths));
        tmpy = nan(1,length(wavelengths));
        % Loop through wavelengths
        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            wavelength = wavelengths(w) * rotor_diameter;
            wavenumber = (2 * pi) / (wavelength);
            amplitude = (wavelength * steepness) / (2 * pi);

            %%% Different geometric scalings
            % Phase difference between front/rear buoys
            buoy_phase = (floater_spacing * cosd(30)) / wavelength;

            % Max geometric pitch angle
            geometric_pitch = rad2deg((2 * amplitude / floater_spacing) * sin((wavenumber * floater_spacing) / 2));


            % Get eigenvectors for one case
            V = decomposition.(farm_spacing).(wave).eigenvectors;
            
            % Mode 1 contributions
            mode1_contrib = V(:,1).^2;
            
            % Quantity 2: purity of mode 1
            mode1_purity = max(mode1_contrib);
            
            % Quantity 3: dominant DOF index of mode 1
            [~, mode1_dof_idx] = max(mode1_contrib);
            
            % If you have names
            dof_names = {'surge','heave','sway','roll','pitch','yaw'};
            mode1_dominant_dof = dof_names{mode1_dof_idx};

            % Plot
            scatter(buoy_phase, mode1_purity, ...
                    50, 'filled', 'markerfacecolor', wave_colors(w,:),...
                    'HandleVisibility', 'off')

            tmpx(w) = buoy_phase;
            tmpy(w) = mode1_purity;
        end

        % Plot connecting line
        [sorted_tmpx, idx] = sort(tmpx, 'ascend');
        P = plot(sorted_tmpx, tmpy(idx), 'color', 'k', ...
                 'linestyle', steepness_linestyle{st}, ...
                 'linewidth', 1.5, ...
                 'HandleVisibility', 'off');
        uistack(P, 'bottom')
    end
    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
end

% Legends
hold on
for w = 1:length(wavelengths)
    plot(nan,nan,'LineWidth',2,'Color',wave_colors(w,:), ...
        'DisplayName',sprintf('$\\lambda=%dD$', wavelengths(w)));
end
plot(nan,nan,'Color','w','DisplayName',' ')
for st = 1:length(wave_steepnesses)
    plot(nan,nan,'k','LineWidth',2,'LineStyle',steepness_linestyle{st}, ...
        'DisplayName',sprintf('$ak=%0.2f$', wave_steepnesses(st)));
end
hold off
leg = legend('Interpreter','latex','Box','off','location', 'northeast'); 
% leg.Layout.Tile = 'east';


linkaxes(h, 'xy')
ylim([0,1])
xlim([0,0.5])
ylabel(tile, 'Maximum Composition: Mode 1', 'interpreter', 'latex')
xlabel(tile, 'Front/Rear Buoy Phase Misalingment: $s_x \mathbin{/} \lambda$', 'interpreter', 'latex')
clc;




%% What DOF contributes the most in Mode 1


clc; close all
figure('color', 'white')
tile = tiledlayout(1,1);

% Loop through farm spacings
% for s = 1:length(farm_spacings)
for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    h(s) = nexttile;
    hold on
    % Loop through wave steepnesses
    % for st = 1:length(wave_steepnesses)
    for st = 1:3
    
        tmpx = nan(1,length(wavelengths));
        tmpy = nan(1,length(wavelengths));
        % Loop through wavelengths
        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            wavelength = wavelengths(w) * rotor_diameter;
            wavenumber = (2 * pi) / (wavelength);
            amplitude = (wavelength * steepness) / (2 * pi);

            %%% Different geometric scalings
            % Phase difference between front/rear buoys
            buoy_phase = (floater_spacing * cosd(30)) / wavelength;

            % Max geometric pitch angle
            geometric_pitch = rad2deg((2 * amplitude / floater_spacing) * sin((wavenumber * floater_spacing) / 2));


            % Get eigenvectors for one case
            V = decomposition.(farm_spacing).(wave).eigenvectors;
            
            % Mode 1 contributions
            mode1_contrib = V(:,1).^2;
            
            % Quantity 2: purity of mode 1
            mode1_purity = max(mode1_contrib);
            
            % Quantity 3: dominant DOF index of mode 1
            [~, mode1_dof_idx] = max(mode1_contrib);
            
            % If you have names
            dof_names = {'surge','heave','sway','roll','pitch','yaw'};
            mode1_dominant_dof = dof_names{mode1_dof_idx};

            % Plot
            scatter(wavelength, mode1_dof_idx, ...
                    50, 'filled', 'markerfacecolor', wave_colors(w,:),...
                    'HandleVisibility', 'off')

            tmpx(w) = wavelength;
            tmpy(w) = mode1_dof_idx;
        end

        % Plot connecting line
        % [sorted_tmpx, idx] = sort(tmpx, 'ascend');
        % P = plot(sorted_tmpx, tmpy(idx), 'color', 'k', ...
        %          'linestyle', steepness_linestyle{st}, ...
        %          'linewidth', 1.5, ...
        %          'HandleVisibility', 'off');
        % uistack(P, 'bottom')
    end
    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
end

% Legends
hold on
for w = 1:length(wavelengths)
    plot(nan,nan,'LineWidth',2,'Color',wave_colors(w,:), ...
        'DisplayName',sprintf('$\\lambda=%dD$', wavelengths(w)));
end
plot(nan,nan,'Color','w','DisplayName',' ')
% for st = 1:length(wave_steepnesses)
%     plot(nan,nan,'k','LineWidth',2,'LineStyle',steepness_linestyle{st}, ...
%         'DisplayName',sprintf('$ak=%0.2f$', wave_steepnesses(st)));
% end
hold off
leg = legend('Interpreter','latex','Box','off','location', 'northeast'); 
% leg.Layout.Tile = 'east';


linkaxes(h, 'xy')
yticks(1:6)
yticklabels(dof_names)
ylim([0,7])
ylabel('Dominant DOF in Mode 1')
ylabel(tile, 'Maximum Composition: Mode 1', 'interpreter', 'latex')
xlabel(tile, 'Front/Rear Buoy Phase Misalingment: $s_x \mathbin{/} \lambda$', 'interpreter', 'latex')
clc;






%% Effective number of active modes

% clc; close all
% figure('color', 'white')
% tile = tiledlayout(1,5);
% 
% % Loop through farm spacings
% for s = 1:length(farm_spacings)
% % for s = 1
%     farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];
% 
%     h(s) = nexttile;
%     hold on
%     % Loop through wave steepnesses
%     % for st = 1:length(wave_steepnesses)
%     for st = 1:3
% 
%         tmpx = nan(1,length(wavelengths));
%         tmpy = nan(1,length(wavelengths));
%         % Loop through wavelengths
%         for w = 1:length(wavelengths)
%             steepness = wave_steepnesses(st);
%             steep = compose('%02d', round(100 * steepness));
%             wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];
% 
%             wavelength = wavelengths(w) * rotor_diameter;
%             wavenumber = (2 * pi) / (wavelength);
%             amplitude = (wavelength * steepness) / (2 * pi);
% 
%             %%% Different geometric scalings
%             % Phase difference between front/rear buoys
%             buoy_phase = rad2deg((2 * pi * floater_spacing) / (wavelength));
% 
%             % Max geometric pitch angle
%             geometric_pitch = rad2deg((2 * amplitude / floater_spacing) * sin((wavenumber * floater_spacing) / 2));
% 
%             lambda = decomposition.(farm_spacing).(wave).eigenvalues;
%             lambda = abs(lambda);
%             Neff = (sum(lambda)^2) / sum(lambda.^2);
% 
%             scatter(buoy_phase, Neff, ...
%                     20, 'filled', 'markerfacecolor', wave_colors(w,:),...
%                     'HandleVisibility', 'off')
% 
%             tmpx(w) = buoy_phase;
%             tmpy(w) = Neff;
%         end
%         [sorted_tmpx, idx] = sort(tmpx, 'ascend');
%         P = plot(sorted_tmpx, tmpy(idx), 'color', 'k', 'linestyle', steepness_linestyle{st}, 'HandleVisibility', 'off');
%         uistack(P, 'bottom')
%     end
%     hold off
%     title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
% end
% 
% % Legends
% hold on
% for w = 1:length(wavelengths)
%     plot(nan,nan,'LineWidth',2,'Color',wave_colors(w,:), ...
%         'DisplayName',sprintf('$\\lambda=%dD$', wavelengths(w)));
% end
% plot(nan,nan,'Color','w','DisplayName',' ')
% for st = 1:length(wave_steepnesses)
%     plot(nan,nan,'k','LineWidth',2,'LineStyle',steepness_linestyle{st}, ...
%         'DisplayName',sprintf('$ak=%0.2f$', wave_steepnesses(st)));
% end
% hold off
% leg = legend('Interpreter','latex','Box','off'); leg.Layout.Tile = 'east';
% 
% 
% linkaxes(h, 'xy')
% ylim([0,6])
% % xlim([0,1])
% ylabel(tile, '\# Active DOF', 'interpreter', 'latex')
% xlabel(tile, '$mod(\lambda \mathbin{/} L_{x}, 1)$', 'interpreter', 'latex')



%% Pitch contribution in mode 1


mode = 1;

clc; close all
figure('color', 'white')
tile = tiledlayout(1,1);

% Loop through farm spacings
% for s = 1:length(farm_spacings)
for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    h(s) = nexttile;
    hold on
    % Loop through wave steepnesses
    % for st = 1:length(wave_steepnesses)
    for st = 1:3
    
        tmpx = nan(1,length(wavelengths));
        tmpy = nan(1,length(wavelengths));
        % Loop through wavelengths
        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            wavelength = wavelengths(w) * rotor_diameter;
            wavenumber = (2 * pi) / (wavelength);
            amplitude = (wavelength * steepness) / (2 * pi);

            %%% Different geometric scalings
            % Phase difference between front/rear buoys
            buoy_phase = floater_spacing / wavelength;

            % Max geometric pitch angle
            geometric_pitch = rad2deg((2 * amplitude / floater_spacing) * sin((wavenumber * floater_spacing) / 2));

            lambda = decomposition.(farm_spacing).(wave).eigenvectors(5, mode).^2;


            scatter(buoy_phase, lambda, ...
                    20, 'filled', 'markerfacecolor', wave_colors(w,:),...
                    'HandleVisibility', 'off')

            % tmpx(w) = mod(wavelength / platform_length_x, 1);
            tmpx(w) = buoy_phase;
            tmpy(w) = lambda;
        end
        [sorted_tmpx, idx] = sort(tmpx, 'ascend');
        P = plot(sorted_tmpx, tmpy(idx), 'color', 'k', 'linestyle', steepness_linestyle{st}, 'HandleVisibility', 'off');
        uistack(P, 'bottom')
    end
    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
end

% Legends
hold on
for w = 1:length(wavelengths)
    plot(nan,nan,'LineWidth',2,'Color',wave_colors(w,:), ...
        'DisplayName',sprintf('$\\lambda=%dD$', wavelengths(w)));
end
plot(nan,nan,'Color','w','DisplayName',' ')
for st = 1:length(wave_steepnesses)
    plot(nan,nan,'k','LineWidth',2,'LineStyle',steepness_linestyle{st}, ...
        'DisplayName',sprintf('$ak=%0.2f$', wave_steepnesses(st)));
end
hold off
leg = legend('Interpreter','latex','Box','off'); leg.Layout.Tile = 'east';


linkaxes(h, 'xy')
ylim([0,1])
xlim([0,0.5])
ylabel(tile,'$v_{\mathrm{pitch},1}^2$','Interpreter','latex')
xlabel(tile, 'Front/Rear Buoy Phase Misalingment: $s \mathbin{/} \lambda$', 'interpreter', 'latex')



%% Largest pitch contribution among first few modes

modes_to_search = 1:3;   % search only energetic modes
pitch_row = 5;           % surge, heave, sway, roll, pitch, yaw

clc; close all
figure('color', 'white')
tile = tiledlayout(1,1);

for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    h(s) = nexttile;
    hold on

    for st = 1:3

        tmpx = nan(1,length(wavelengths));
        tmpy = nan(1,length(wavelengths));

        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            wavelength = wavelengths(w) * rotor_diameter;
            wavenumber = (2 * pi) / wavelength;
            amplitude = (wavelength * steepness) / (2 * pi);

            %%% Different geometric scalings
            % Phase difference proxy
            buoy_phase = floater_spacing / wavelength;

            % Max geometric pitch angle
            geometric_pitch = rad2deg((2 * amplitude / floater_spacing) * ...
                                      sin((wavenumber * floater_spacing) / 2));

            % Get eigenvectors
            V = decomposition.(farm_spacing).(wave).eigenvectors;

            % Search only first few modes
            pitch_contribs = V(pitch_row, modes_to_search).^2;

            % Find maximum pitch contribution and corresponding mode
            [max_pitch_contrib, local_idx] = max(pitch_contribs);
            best_mode = modes_to_search(local_idx);

            % Find maximum pitch wighted by modal contribution
            % eigvals = decomposition.(farm_spacing).(wave).fractional_energy;
            % pitch_weighted = eigvals(modes_to_search)' .* (V(pitch_row, modes_to_search).^2);
            % [max_pitch_contrib, local_idx] = max(pitch_weighted);
            % best_mode = modes_to_search(local_idx);

            

            % Plot point
            scatter(buoy_phase, max_pitch_contrib, ...
                    50, 'filled', ...
                    'markerfacecolor', wave_colors(w,:), ...
                    'HandleVisibility', 'off')

            % Annotate with mode number
            % text(buoy_phase + 0.005, max_pitch_contrib, ...
            %      sprintf('m%d', best_mode), ...
            %      'FontSize', 8, ...
            %      'Color', 'k', ...
            %      'HorizontalAlignment', 'left', ...
            %      'VerticalAlignment', 'middle');

            tmpx(w) = buoy_phase;
            tmpy(w) = max_pitch_contrib;
        end

        [sorted_tmpx, idx] = sort(tmpx, 'ascend');
        P = plot(sorted_tmpx, tmpy(idx), ...
                 'color', 'k', ...
                 'linestyle', steepness_linestyle{st}, ...
                 'linewidth', 1.5, ...
                 'HandleVisibility', 'off');
        uistack(P, 'bottom')
    end

    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
end

% Legends
hold on
for w = 1:length(wavelengths)
    plot(nan,nan,'LineWidth',2,'Color',wave_colors(w,:), ...
        'DisplayName',sprintf('$\\lambda=%dD$', wavelengths(w)));
end
plot(nan,nan,'Color','w','DisplayName',' ')
for st = 1:length(wave_steepnesses)
    plot(nan,nan,'k','LineWidth',2,'LineStyle',steepness_linestyle{st}, ...
        'DisplayName',sprintf('$ak=%0.2f$', wave_steepnesses(st)));
end
hold off
leg = legend('Interpreter','latex','Box','off');
leg.Layout.Tile = 'east';

linkaxes(h, 'xy')
ylim([0,1])
xlim([0,0.5])
ylabel(tile,'Largest pitch contribution in searched modes','Interpreter','latex')
xlabel(tile,'Front/Rear Buoy Phase Misalignment: $s/\lambda$','Interpreter','latex')


%% Distribution within mode 1 (squared loadings)

mode = 1;
clc; close all
figure('color','white')
tile = tiledlayout(1,5);
sgtitle(sprintf('Mode %d composition', mode))

for s = 1:length(farm_spacings)
    farm_spacing = ['SX', num2str(farm_spacings(s)*10)];
    h(s) = nexttile; hold on

    % for st = 1:length(wave_steepnesses)
    for st = 3
        for w = 1:length(wavelengths)

            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100*steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            v = decomposition.(farm_spacing).(wave).eigenvectors(:,mode);
            contrib = (v.^2) / sum(v.^2);

            plot(1:6, contrib, 'LineWidth', 1.5, ...
                 'Color', wave_colors(w,:), ...
                 'LineStyle', steepness_linestyle{st}, ...
                 'HandleVisibility','off');
        end
    end

    xticks(1:6); 
    xticklabels(dof_names);
    ylim([0 1]); 
    xlim([0.5 6.5])
    % grid on; box on
    title(sprintf('$S_x=%1.1fD$', farm_spacings(s)), 'Interpreter','latex')
end

% Legends
hold on
for w = 1:length(wavelengths)
    plot(nan,nan,'LineWidth',2,'Color',wave_colors(w,:), ...
        'DisplayName',sprintf('$\\lambda=%dD$', wavelengths(w)));
end
plot(nan,nan,'Color','w','DisplayName',' ')
for st = 1:length(wave_steepnesses)
    plot(nan,nan,'k','LineWidth',2,'LineStyle',steepness_linestyle{st}, ...
        'DisplayName',sprintf('$ak=%0.2f$', wave_steepnesses(st)));
end
hold off
leg = legend('Interpreter','latex','Box','off'); leg.Layout.Tile = 'east';

ylabel(tile,'Contribution','Interpreter','latex')



%% Test spider plot


mode = 1;
clc; close all
figure('color','white')
tile = tiledlayout(1,1);
sgtitle(sprintf('Mode %d composition', mode))

% for s = 1:length(farm_spacings)
for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s)*10)];
    h(s) = nexttile; hold on

    % for st = 1:length(wave_steepnesses)
    tmp = nan(4,6);
    for st = 3
        for w = 1:length(wavelengths)

            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100*steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            v = decomposition.(farm_spacing).(wave).eigenvectors(:,mode);
            contrib = (v.^2) / sum(v.^2);

            tmp(w,:) = contrib.';
        end
    end

    spider_plot_R2019b(tmp, ...
            'AxesInterval', 1,...
            'AxesPrecision', 0, ...
            'AxesLimits', [0,0,0,0,0,0;1,1,1,1,1,1], ...
            'AxesLabels', {'Surge: $x$', 'Heave: $y$', 'Sway: $z$', 'Roll: $\phi$', 'Pitch: $\theta$', 'Yaw: $\psi$'}, ...
            'AxesInterpreter', 'latex', ...
            'AxesLabelsEdge', 'none', ...
            'LineWidth', 2 * ones(1,4), ...
            'AxesDisplay', 'one',...
            'FillOption', 'on', ...
            'FillTransparency', 0.2,...
            'Color', wave_colors, ...
            'AxesDisplay', 'one', ...
            'AxesOffset', 1,...
            'LabelFontSize', 8);

    title(sprintf('$S_x=%1.1fD, ak = %1.2f$', farm_spacings(s), steepness), 'Interpreter','latex')
end

% Legends
hold on
for w = 1:length(wavelengths)
    plot(nan,nan,'LineWidth',2,'Color',wave_colors(w,:), ...
        'DisplayName',sprintf('$\\lambda=%dD$', wavelengths(w)));
end
hold off
leg = legend('Interpreter','latex','Box','off', 'orientation', 'horizontal'); leg.Layout.Tile = 'south';



