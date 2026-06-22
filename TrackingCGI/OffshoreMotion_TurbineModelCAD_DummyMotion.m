
%% Motion movie test (Better saving)
% Saves each animation frame as a PNG, assembles the PNGs into an MP4
% using MATLAB VideoWriter, and then deletes the temporary frame folder.

clc
close all

% Paths
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

% Constants
tickFontSize = 8;
labelFontSize = 10;

markerSize = 15;
turbineColor = [0.80 0.80 0.10];

% Number of recent frames retained in the trail
trailLength = 30;
shadowColormap = 'BuPu';

% Frames per second
frameRate = 30;   

% Seconds
duration  = 5;  

% Motion frequency
fMotion = 0.4;    % Hz

% Translational amplitudes [mm]
A_surge = 30;
A_sway  = 5;
A_heave = 4;

% Rotational amplitudes [deg]
A_roll  = 3;
A_pitch = 5;
A_yaw   = 10;

% Saves
videoSavePath = '/Users/zeinsadek/Desktop/Experiments/Offshore/Tracking/Processing/TrackingCGI';
videoSaveName = 'OffshoreMotion_DummyMotion.mp4';

% PNG export resolution
frameResolution = 300;


%% Read and prepare STL model
modelData = stlread('/Users/zeinsadek/Desktop/Experiments/Offshore/Tracking/Processing/TrackingCGI/FOWT.stl');

%%% ORIGIN
% From arbitrary assembly origin to one corner of the assembly
shiftx = -58.110;
shifty =  52.426;
shiftz = -149.767;

% Additional shifts to place the origin at the tower center
subtleShiftx = -9.4/2;
subtleShifty =  11.4/2;

x_origin = -shiftx - subtleShiftx;
y_origin = -shifty - subtleShifty;
z_origin = -shiftz;

% Place origin at tower base
Vfaces = modelData.ConnectivityList;
Vcad   = modelData.Points;

Vcentered = Vcad - [x_origin, y_origin, z_origin];

% Correct coordinate alignment
Vdata = zeros(size(Vcentered));
Vdata(:,1) = -Vcentered(:,2);
Vdata(:,2) =  Vcentered(:,1);
Vdata(:,3) =  Vcentered(:,3);

% Create figure
fig = figure( ...
    'Color', 'white', ...
    'Units', 'centimeters', ...
    'Position', [10, 10, 10, 10]);

ax = axes('Parent', fig);
hold(ax, 'on')

fig.Renderer = 'opengl';
ax.SortMethod = 'depth';
set(ax, 'FontSize', tickFontSize)
set(ax, 'TickLabelInterpreter', 'latex');

% Reference origin
scatter3(ax, 0, 0, 0, ...
    markerSize, ...
    'filled', ...
    'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'none');

% Transformation container
motionTransform = hgtransform('Parent', ax);

% STL model
hModel = patch( ...
    'Faces', Vfaces, ...
    'Vertices', Vdata, ...
    'FaceColor', turbineColor, ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'flat', ...
    'Parent', motionTransform);

axis(ax, 'equal')
view(ax, -45, 30)
grid(ax, 'on')

xlabel(ax, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(ax, '$z$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
zlabel(ax, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)

camlight(ax, 'headlight')
material(ax, 'dull')

% Reverse MATLAB Y to match the experimental spanwise convention
set(ax, 'YDir', 'reverse')




%%% Motion-shadow points
% Each row is one reference point attached to the undeformed turbine.
% Coordinates must be in the same coordinate system as Vdata.
shadowPoints0 = [ ...
      17.5,   0.0, 100.0;   % Hub
    -41.9,  72.5,   0.0;   % Front-left buoy
    -41.9, -72.5,   0.0;   % Front-right buoy
     83.7,   0.0,   0.0];  % Rear buoy

nShadowPoints = size(shadowPoints0, 1);

% Dimensions:
% trailLength x nShadowPoints x 3
pointHistory = nan(trailLength, nShadowPoints, 3);

% Colormap for oldest-to-newest points
trailColors = slanCM(shadowColormap, trailLength);

% One scatter object for all trail points
hShadow = scatter3(ax, ...
    nan, nan, nan, ...
    markerSize, ...
    nan, ...
    'filled', ...
    'MarkerEdgeColor', 'none');

colormap(ax, trailColors)




%%% Dummy six-DOF motion
nFrames   = frameRate * duration;

% Exclude duplicated endpoint
t = (0:nFrames-1) / frameRate;

% Dummy time series
surgeData = A_surge * sin(2*pi*fMotion*t);
swayData  = A_sway  * sin(2*pi*fMotion*t + pi/3);
heaveData = A_heave * sin(2*pi*fMotion*t + pi/2);

rollData  = A_roll  * sin(2*pi*fMotion*t + pi/4);
pitchData = A_pitch * sin(2*pi*fMotion*t);
yawData   = A_yaw   * sin(2*pi*fMotion*t + 2*pi/3);

% Fix axis limits

modelExtent = max(Vdata, [], 1) - min(Vdata, [], 1);
margin = 0.15 * max(modelExtent);

xlim(ax, [ ...
    min(Vdata(:,1)) - A_surge - margin, ...
    max(Vdata(:,1)) + A_surge + margin]);

ylim(ax, [ ...
    min(Vdata(:,2)) - A_sway - margin, ...
    max(Vdata(:,2)) + A_sway + margin]);

zlim(ax, [ ...
    min(Vdata(:,3)) - A_heave - margin, ...
    max(Vdata(:,3)) + A_heave + margin]);

axis(ax, 'manual')

% Time indicator
hTime = title(ax, ...
    sprintf('$t = %.2f$ s', t(1)), ...
    'Interpreter', 'latex');

% Prepare temporary frame folder
videoFile = fullfile(videoSavePath, videoSaveName);
tmpFolder = fullfile(videoSavePath, 'FOWT_temp_frames');

% Remove old temporary folder if it exists
if exist(tmpFolder, 'dir')
    rmdir(tmpFolder, 's');
end
mkdir(tmpFolder)

% Zero-padded filenames ensure correct chronological sorting
nDigits = max(4, ceil(log10(nFrames + 1)));

fprintf('Saving temporary frames to:\n%s\n\n', tmpFolder);

% Transformation origin
rotationOrigin = [0, 0, 0];

T_toOrigin = eye(4);
T_toOrigin(1:3,4) = -rotationOrigin(:);

T_fromOrigin = eye(4);
T_fromOrigin(1:3,4) = rotationOrigin(:);

% Animation loop

for n = 1:nFrames

    % Translations [mm]

    % MATLAB X = experimental x
    surge = surgeData(n);

    % MATLAB Y = experimental z
    sway = -swayData(n);

    % MATLAB Z = experimental y
    heave = heaveData(n);

    % Rotations [rad]

    % Roll about MATLAB X
    roll = -deg2rad(rollData(n));

    % Pitch about MATLAB Y
    pitch = -deg2rad(pitchData(n));

    % Yaw about MATLAB Z
    yaw = -deg2rad(yawData(n));

    % Rotation matrices

    Rx = [ ...
        1,          0,           0, 0;
        0,  cos(roll), -sin(roll), 0;
        0,  sin(roll),  cos(roll), 0;
        0,          0,           0, 1];

    Ry = [ ...
         cos(pitch), 0, sin(pitch), 0;
                  0, 1,          0, 0;
        -sin(pitch), 0, cos(pitch), 0;
                  0, 0,          0, 1];

    Rz = [ ...
        cos(yaw), -sin(yaw), 0, 0;
        sin(yaw),  cos(yaw), 0, 0;
               0,         0, 1, 0;
               0,         0, 0, 1];

    % Apply roll, then pitch, then yaw
    R = Rz * Ry * Rx;

    % Translation matrix
    T_motion = eye(4);
    T_motion(1:3,4) = [surge; sway; heave];

    % Full rigid-body transformation
    M = T_motion * T_fromOrigin * R * T_toOrigin;
    motionTransform.Matrix = M;

    % Transform shadow points
    P0_h = [shadowPoints0, ones(nShadowPoints,1)];
    P_h = (M * P0_h.').';
    P   = P_h(:,1:3);

    % Shift history backward and append newest positions
    pointHistory(1:end-1,:,:) = pointHistory(2:end,:,:);
    pointHistory(end,:,:) = reshape(P, 1, nShadowPoints, 3);

    % Keep valid rows during startup
    validTime = ~isnan(pointHistory(:,1,1));
    historyValid = pointHistory(validTime,:,:);
    nValid = size(historyValid,1);

    % Flatten as:
    % time 1, all tracked points
    % time 2, all tracked points
    % ...
    Pplot = reshape(permute(historyValid, [2, 1, 3]), [], 3);

    % Use one unique color per retained time step
    Cvalid = trailColors(end-nValid+1:end, :);

    % Repeat each age color for all tracked locations
    Cplot = repelem(Cvalid, nShadowPoints, 1);

    % Update scatter object
    hShadow.XData = Pplot(:,1);
    hShadow.YData = Pplot(:,2);
    hShadow.ZData = Pplot(:,3);
    hShadow.CData = Cplot;

    % Update time and render
    hTime.String = sprintf('$t = %.2f$ s', t(n));

    drawnow

    % Save current frame as PNG
    frameFile = fullfile(tmpFolder, ...
        sprintf(['frame_%0', num2str(nDigits), 'd.png'], n));

    exportgraphics(fig, frameFile, ...
        'Resolution', frameResolution, ...
        'BackgroundColor', 'white');

    fprintf('Saved frame %d of %d\n', n, nFrames);
end




%%% Assemble PNG frames into MP4
fprintf('\nAssembling video...\n');

% Get list of frame images
frames = dir(fullfile(tmpFolder, '*.png'));

% Sort frames by name if their names are zero-padded
[~, idx] = sort({frames.name});
frames = frames(idx);

% Read first image to establish the video dimensions
firstFrame = imread(fullfile(tmpFolder, frames(1).name));

originalHeight = size(firstFrame, 1);
originalWidth  = size(firstFrame, 2);

% Force dimensions to be multiples of 8
videoHeight = floor(originalHeight / 8) * 8;
videoWidth  = floor(originalWidth  / 8) * 8;

fprintf('Original first-frame size: %d x %d\n', ...
    originalHeight, originalWidth);

fprintf('Video frame size: %d x %d\n', ...
    videoHeight, videoWidth);

% Create video
v = VideoWriter(videoFile, 'MPEG-4');
v.FrameRate = frameRate;
v.Quality = 100;
open(v);

for k = 1:numel(frames)

    framePath = fullfile(tmpFolder, frames(k).name);
    img = imread(framePath);

    % Handle indexed, grayscale, or RGBA images
    if size(img, 3) == 1
        img = repmat(img, 1, 1, 3);
    elseif size(img, 3) == 4
        img = img(:, :, 1:3);
    end

    currentHeight = size(img, 1);
    currentWidth  = size(img, 2);

    fprintf('Frame %d original size: %d x %d\n', ...
        k, currentHeight, currentWidth);

    % Resize every frame to the exact same dimensions
    if currentHeight ~= videoHeight || currentWidth ~= videoWidth
        img = imresize(img, [videoHeight, videoWidth]);
    end

    % Safety check
    assert(size(img,1) == videoHeight, ...
        'Incorrect frame height at frame %d.', k);

    assert(size(img,2) == videoWidth, ...
        'Incorrect frame width at frame %d.', k);

    writeVideo(v, img);

    fprintf('Added frame %d of %d\n', k, numel(frames));
end

close(v);
fprintf('Video saved successfully:\n%s\n', videoFile);

% Clean up temporary images
disp('Movie generation complete.')
fprintf('Saved video: %s\n', videoFile);

if exist(tmpFolder, 'dir')
    rmdir(tmpFolder, 's');
end
fprintf('Temporary frames folder deleted!\n');
close all