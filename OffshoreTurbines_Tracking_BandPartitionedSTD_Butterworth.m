%% Looking at band-limited RMS/STD of motion
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
% COMPUTE BAND-LIMITED RMS OF EACH DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = 30;
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
        % Get wave and wave frequency
        wave = waves{w};
        split_wave = split(wave, '_');
        wavelength_name = split_wave{1};

        if ~strcmp(wavelength_name, 'LM0')
            fw = forcing_frequencies.(wavelength_name);

            % record length based fmin (or just use 1/T)
            tvec = tracking.(farm_spacing).(wave)(t).time;
            Trec = tvec(end) - tvec(1);
            fmin = max(1/Trec, 0.01);   % keep your 0.01 floor if you want
            % fmin = 1/Trec;
            
            % Contiguous bands (no gaps)
            bands.LF   = [fmin, 0.5*fw];
            bands.Wave = [0.5*fw, 1.5*fw];
            bands.HF   = [1.5*fw, 0.49*Fs];

    
            % Loop through DOFs
            for d = 1:length(DOFs)
                DOF = DOFs{d};
    
                % Loop through turbines
                for t = 1:length(turbines)
                    % Signal
                    q = tracking.(farm_spacing).(wave)(t).(DOF);
    
                    % Compute low frequency portion
                    [rLF, qLF]     = bandlimited_rms(q, Fs, bands.LF(1), bands.LF(2), 'constant');
                    
                    % Compute wave frequency portion
                    [rWave, qWave] = bandlimited_rms(q, Fs, bands.Wave(1), bands.Wave(2), 'constant');
    
                    % Compute high frequency portion
                    [rHF, qHF]     = bandlimited_rms(q, Fs, bands.HF(1), bands.HF(2), 'constant');
    
    
                    % Save
                    bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).LF = rLF;
                    bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).WF = rWave;
                    bandfilteredDeviations.(farm_spacing).(wave)(t).(DOF).HF = rHF;
    
                end
            end
        end
    end
end



%% CHeck that sum of square of bands is close to the full RMS

clc;
spacing = 'SX50';
wave = 'LM5_AK12';
turbine = 8;
DOF = 'pitch_kal';

full = deviations.(spacing).(wave)(turbine).(DOF);
LF = bandfilteredDeviations.(spacing).(wave)(turbine).(DOF).LF;
WF = bandfilteredDeviations.(spacing).(wave)(turbine).(DOF).WF;
HF = bandfilteredDeviations.(spacing).(wave)(turbine).(DOF).HF;
partition = [LF, WF, HF];

% Sum of squares of partitioned signals
ss = sum(partition.^2);
total_squared = full^2;


disp(partition)
disp(full)
disp(ss)
disp(ss / total_squared)

%% testing












%% Functions

function [rmsBand, qBand] = bandlimited_rms(q, Fs, f1, f2, detrendMode)
% bandlimited_rms  Compute RMS of q in the frequency band [f1,f2] (Hz)
%
% Inputs:
%   q           [Nx1] signal
%   Fs          sampling frequency [Hz]
%   f1,f2       band edges [Hz] (0 < f1 < f2 < Fs/2)
%   detrendMode 'constant' (mean removal) or 'linear'
%
% Outputs:
%   rmsBand     RMS of bandpassed signal
%   qBand       bandpassed signal (zero-phase)

    arguments
        q (:,1) double
        Fs (1,1) double {mustBePositive}
        f1 (1,1) double {mustBeNonnegative}
        f2 (1,1) double {mustBePositive}
        detrendMode (1,:) char = 'constant'
    end

    if f2 >= Fs/2
        error('f2 must be < Fs/2. Got f2=%.3f, Fs/2=%.3f', f2, Fs/2);
    end
    if f1 <= 0
        error('Use a highpass cutoff > 0 for bandpass RMS. Got f1=%.3f', f1);
    end

    q0 = detrend(q, detrendMode);

    % 4th-order Butterworth bandpass
    [b,a] = butter(4, [f1 f2]/(Fs/2), 'bandpass');
    qBand = filtfilt(b,a,q0);   % zero-phase (no lag)

    rmsBand = sqrt(mean(qBand.^2));
end





