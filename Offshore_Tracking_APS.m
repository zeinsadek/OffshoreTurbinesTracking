%%% Offshore Tracking Signals

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/Offshore_Functions/');
fprintf('All Paths Imported...\n\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
path = '/Users/zeinsadek/Desktop/APS2024_Data/Zein_APS_Tracking/Tracking/FWF_WT60_SX50_AG0.mat';
out_path = '/Users/zeinsadek/Desktop/APS2024_Data/Tracking_Results';
data = load(path);
data = data.('WT60_SX50_AG0'); 
center_turbines = [2,5,8,11];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECT WAVE CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wave = 'LM5_AK12';
signals = data.(wave);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT TURBINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(center_turbines)
    turbine = strcat('turbine_', num2str(center_turbines(i)));
    center_signals.(turbine) = signals.(turbine);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sliding mean window size
window = 6;

% Colors
colors = {'#710627', '#EF6461', '#9CE37D', '#8ACDEA'};

% Plot details
linewidth = 1.5;
marker_size = 15;

% Plot
ax = figure();
tiledlayout(4,1)
sgtitle(wave, 'interpreter', 'none')
for i = 1:length(center_turbines)
    nexttile
    turbine = strcat('turbine_', num2str(center_turbines(i)));

    hold on
    plot(center_signals.(turbine).time, movmean(center_signals.(turbine).pitch, window), ...
         'Color', colors{i}, 'LineWidth', linewidth)
    
    scatter(center_signals.(turbine).time, movmean(center_signals.(turbine).pitch, window), ...
            marker_size, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors{i})
    hold off

    title(strcat("Row", " ", num2str(i)), 'Interpreter', 'None')
    xlim([10, 20])
    ylim([-5, 20])

    if i == 2
        ylabel('Pitch [Degrees]')
    end

    if i == 4
        xlabel('Time [Seconds]')
    end
end

% Save fig
% exportgraphics(ax, fullfile(out_path, strcat(wave, '_Pitch_Sample.png')), 'Resolution', 200)


%% Plot w signals overlayed

% Sliding mean window size
window = 6;

% Colors
colors = {'#710627', '#EF6461', '#9CE37D', '#8ACDEA'};

% Plot details
linewidth = 2.5;
marker_size = 30;

% Plot
ax = figure();
sgtitle(wave, 'interpreter', 'none')
hold on
for i = 1:length(center_turbines)
    turbine = strcat('turbine_', num2str(center_turbines(i)));
    plot(center_signals.(turbine).time, movmean(center_signals.(turbine).pitch, window), ...
         'Color', colors{i}, 'LineWidth', linewidth, 'DisplayName', strcat('Row', " ", num2str(i)))
    
    scatter(center_signals.(turbine).time, movmean(center_signals.(turbine).pitch, window), ...
            marker_size, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors{i}, ...
            'HandleVisibility', 'off') 
end
hold off
xlim([16, 17])  
ylim([-10, 25])
legend()
ylabel('Pitch [Degrees]')
xlabel('Time [Seconds]')

% Save fig
% exportgraphics(ax, fullfile(out_path, strcat(wave, '_Pitch_Overlayed.png')), 'Resolution', 200)

%% FFT

% LM33
if strcmp(wave, 'LM33_AK12') == 1
    wave_freq = 1.7651;
end

% LM50
if strcmp(wave, 'LM5_AK12') == 1
    wave_freq = 1.4273;
end

ax = figure();
hold on
for i = 1:length(center_turbines)

    turbine = strcat('turbine_', num2str(center_turbines(i)));
    time    = center_signals.(turbine).time;
    S       = movmean(center_signals.(turbine).pitch, 5);
    S       = fillmissing(S, 'previous');
    S       = S - mean(S);
    
    Fs = 1/mean(diff(time), 'omitnan');                         
    T  = 1/Fs;           
    L  = length(S);           
    
    Y  = fft(S);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs/L*(0:(L/2));
    
    plot(f, P1, "LineWidth", 3, "Color", colors{i}) 

end
hold off
xline(0.5 * wave_freq, 'color', 'black', 'linestyle', '--')
xline(wave_freq, 'color', 'black', 'linestyle', '--')
xline(1.5 * wave_freq, 'color', 'black', 'linestyle', '--')
xline(2 * wave_freq, 'color', 'black', 'linestyle', '--')
title(wave, 'Interpreter', 'none')
xlabel("f [Hz]")
ylabel("|Pitch| [Degrees]")
xlim([0,5])
exportgraphics(ax, fullfile(out_path, strcat(wave, '_FFT.png')), 'Resolution', 200)

%% Cross Correlation

window = 5;
ax = figure(Position=[200,200,800,400]);
hold on
for i = 1:length(center_turbines) - 1
    signal1 = movmean(center_signals.('turbine_2').pitch, window);
    signal1 = fillmissing(signal1, 'previous');
    signal1 = signal1 - mean(signal1);

    turbine = strcat('turbine_', num2str(center_turbines(i + 1)));
    signal2 = movmean(center_signals.(turbine).pitch, window);
    signal2 = fillmissing(signal2, 'previous');
    signal2 = signal2 - mean(signal2);

    [cor, lags] = xcorr(signal1, signal2, 20, 'normalized');
    lags = (lags * T) / (1/wave_freq);
    plot(lags, cor, 'DisplayName', strcat('Row 1 X Row', " " , num2str(i + 1)), ...
         'linewidth', 3, 'Color', colors{i + 1})

end
hold on
ylim([-1, 1])
xlim([0, 1])
title(wave, 'Interpreter','none')
xlabel('$t / T_{wave}$', 'interpreter', 'latex')
ylabel('Normalized Correlation')
legend()
exportgraphics(ax, fullfile(out_path, strcat(wave, '_Pitch_Cross_Correlation.png')), 'Resolution', 200)









