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

% Platform length (x)
platform_length_x = 0.1 * cosd(30);
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
surge_norm = surge / platform_length_x;
heave_norm = heave / platform_length_x;
sway_norm  = sway  / platform_length_x;

% Normalize rotations by displacement at hub, over platform length
roll_norm  = (deg2rad(roll) * hub_height_from_surface) / platform_length_x;
pitch_norm = (deg2rad(pitch) * hub_height_from_surface) / platform_length_x;
yaw_norm   = (deg2rad(yaw) * hub_height_from_surface) / platform_length_x;


% Assemble data matrix
X = [surge_norm(:), heave_norm(:), sway_norm(:), ...
     roll_norm(:),  pitch_norm(:), yaw_norm(:)];
X = X - mean(X, 1, 'omitnan');


% PCA
% C = cov(X);            
% [V,D] = eig(C);
% 
% lambda = diag(D);
% [lambda,idx] = sort(lambda,'descend');
% V = V(:,idx);
% 
% explained = lambda/sum(lambda);



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
title(sprintf('6-DOF PCA explained variance (%s, %s, turbine %d)', farm_spacing, wave, turbine));
xlim([0.5 6.5]);

% Optional: cumulative line
% hold on;
% plot(1:6, 100*cumsum(explained), '-o', 'LineWidth', 1.5);
% legend('Individual','Cumulative','Location','best');

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


% 3) Heatmap: covariance matrix
% figure('Color','w');
% imagesc(Cc);
% axis equal tight;
% colorbar;
% set(gca,'XTick',1:6,'XTickLabel',dof_names,'YTick',1:6,'YTickLabel',dof_names);
% title('Covariance matrix cov(X)');
% % Optional numbers on cells:
% for i=1:6
%     for j=1:6
%         text(j,i,sprintf('%.3g',Cc(i,j)),'HorizontalAlignment','center');
%     end
% end

% 4) Heatmap: correlation matrix
figure('Color','w');
imagesc(corr(X));
axis equal tight;
colorbar;
caxis([-1 1]);
set(gca,'XTick',1:6,'XTickLabel',dof_names,'YTick',1:6,'YTickLabel',dof_names);
title('Correlation matrix corr(X)');
% Optional numbers on cells:
for i=1:6
    for j=1:6
        text(j,i,sprintf('%.2f',Rc(i,j)),'HorizontalAlignment','center');
    end
end

% 5) Modal coordinates (time series): eta_k(t) = X * v_k
% eta = Xc * V;   % N x 6, modal coordinates in PCA basis
% 
% figure('Color','w');
% for k = 1:nModesToShow
%     subplot(nModesToShow,1,k);
%     plot(eta(:,k),'LineWidth',1.0);
%     grid on; box on;
%     ylabel(sprintf('\\eta_%d',k));
%     if k==1
%         title('Modal coordinates (projection of motion onto PCA modes)');
%     end
% end
% xlabel('Sample index');

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
title('Within-mode DOF contributions (squared loadings)');







%% Loop through cases and plot useful metrics


turbine = centers(3);
farm_spacing = 'SX50';
wave = 'LM5_AK12';


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
            % % Translations by platform length
            % surge_norm = surge / platform_length_x;
            % heave_norm = heave / platform_length_x;
            % sway_norm  = sway  / platform_length_x;
            % 
            % % Normalize rotations by displacement at hub, over platform length
            % roll_norm  = (deg2rad(roll) * hub_height_from_surface) / platform_length_x;
            % pitch_norm = (deg2rad(pitch) * hub_height_from_surface) / platform_length_x;
            % yaw_norm   = (deg2rad(yaw) * hub_height_from_surface) / platform_length_x;

            % Normalize: Wave
            % Translations by platform length
            surge_norm = surge / amplitude;
            heave_norm = heave / amplitude;
            sway_norm  = sway  / amplitude;
            
            % Normalize rotations by displacement at hub, over platform length
            roll_norm  = deg2rad(roll) / steepness;
            pitch_norm = deg2rad(pitch) / steepness;
            yaw_norm   = deg2rad(yaw) / steepness;
            
            
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
tile = tiledlayout(1,5);

% Loop through farm spacings
for s = 1:length(farm_spacings)
% for s = 1
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

            % plot(1:6, decomposition.(farm_spacing).(wave).fractional_energy, ...
            %      'linewidth', 2, 'color', wave_colors(w,:))

            plot(100 * cumsum(decomposition.(farm_spacing).(wave).fractional_energy), 1:6, ...
                 'linewidth', 2, 'color', wave_colors(w,:), 'linestyle', steepness_linestyle{st}, ...
                 'HandleVisibility', 'off')

        end
    end
    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
    yticks(1:6)
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

leg = legend('box', 'off', 'interpreter', 'latex', 'orientation', 'vertical');
leg.Layout.Tile = 'east';


linkaxes(h, 'xy')
ylabel(tile, 'Modes', 'interpreter', 'latex')
xlabel(tile, 'Cumulative Energy Represented', 'interpreter', 'latex')



%% Look at the energy contained within the first mode vs wave parameters


clc; close all
figure('color', 'white')
tile = tiledlayout(1,5);


% Loop through farm spacings
for s = 1:length(farm_spacings)
% for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    h(s) = nexttile;
    hold on
    % Loop through wave steepnesses
    % for st = 1:length(wave_steepnesses)
    for st = 1:3
    
        tmp = nan(1,length(wavelengths));
        % Loop through wavelengths
        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            wavelength = wavelengths(w) * rotor_diameter;
            amplitude = (wavelength * steepness) / (2 * pi);

            scatter(wavelength / platform_length_x, decomposition.(farm_spacing).(wave).fractional_energy(1), ...
                    20, 'filled', ...
                    'HandleVisibility', 'off')

            tmp(w) = decomposition.(farm_spacing).(wave).fractional_energy(1);
        end
        plot(wavelengths * rotor_diameter / platform_length_x, tmp)
    end
    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
end

linkaxes(h, 'xy')
ylim([0,1])
xlim([3,9])
ylabel(tile, 'P1', 'interpreter', 'latex')
xlabel(tile, '$\lambda \mathbin{/} L_{x}$', 'interpreter', 'latex')


%% Effective number of active modes

clc; close all
figure('color', 'white')
tile = tiledlayout(1,5);

mode = 1;

% Loop through farm spacings
for s = 1:length(farm_spacings)
% for s = 1
    farm_spacing = ['SX', num2str(farm_spacings(s) * 10)];

    h(s) = nexttile;
    hold on
    % Loop through wave steepnesses
    % for st = 1:length(wave_steepnesses)
    for st = 1:3
    
        tmp = nan(1,length(wavelengths));
        % Loop through wavelengths
        for w = 1:length(wavelengths)
            steepness = wave_steepnesses(st);
            steep = compose('%02d', round(100 * steepness));
            wave = ['LM', num2str(wavelengths(w)), '_AK', steep{1}];

            wavelength = wavelengths(w) * rotor_diameter;
            amplitude = (wavelength * steepness) / (2 * pi);

            eigenvectors = decomposition.(farm_spacing).(wave).eigenvectors(:, mode);

            effective_num_modes = ((sum(eigenvectors).^2) / (sum(eigenvectors.^2)));

            scatter(wavelength / platform_length_x, effective_num_modes, ...
                    20, 'filled', ...
                    'HandleVisibility', 'off')

            tmp(w) = effective_num_modes;
        end
        plot(wavelengths * rotor_diameter / platform_length_x, tmp)
    end
    hold off
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
end

linkaxes(h, 'xy')
ylim([0,6])
xlim([3,9])
ylabel(tile, '\# Active DOF', 'interpreter', 'latex')
xlabel(tile, '$\lambda \mathbin{/} L_{x}$', 'interpreter', 'latex')



%% Distributiom within mode 1

mode = 1;

clc; close all
figure('color', 'white')
tile = tiledlayout(1,5);
sgtitle(sprintf('Mode %1.0f', mode))

% Loop through farm spacings
for s = 1:length(farm_spacings)
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

            wavelength = wavelengths(w) * rotor_diameter;
            amplitude = (wavelength * steepness) / (2 * pi);

            eigenvectors = decomposition.(farm_spacing).(wave).eigenvectors(:, mode).^2;
            contrib = eigenvectors/sum(eigenvectors);

            scatter(1:6, contrib, 20, 'filled', ...
                    'markerfacecolor', wave_colors(w,:), ...
                    'HandleVisibility', 'off')

            plot(1:6, contrib, 'color', wave_colors(w,:), 'linewidth', 2, 'HandleVisibility', 'off')
        end
        
    end
    hold off
    % Set the x-axis tick labels
    xticks(1:6); 
    xticklabels(dof_names);
    title(sprintf('$S_{x} = %1.1fD$', farm_spacings(s)), 'interpreter', 'latex')
end

hold on
% Legend for wave colors
for w = 1:length(wavelengths)
    label = sprintf('$\\lambda = %sD$', num2str(wavelengths(w)));
    plot(nan, nan, 'linewidth', 2, 'color', wave_colors(w,:), 'displayname', label)
end
hold off

leg = legend('box', 'off', 'interpreter', 'latex', 'orientation', 'vertical');
leg.Layout.Tile = 'east';

linkaxes(h, 'xy')
xlim([0.5, 6.5])
ylabel(tile, 'Percentage Contribution', 'interpreter', 'latex')
% xlabel(tile, '$\lambda \mathbin{/} L_{x}$', 'interpreter', 'latex')








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

