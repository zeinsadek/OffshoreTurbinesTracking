%% Motion movie test (Better saving)
% Saves each animation frame as a PNG, assembles the PNGs into an MP4
% using MATLAB VideoWriter, and then deletes the temporary frame folder.

clc; close all; clear

% Paths
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
offshore_path = "/Users/zeinsadek/Desktop/Experiments/Offshore";
tracking_path = fullfile(offshore_path, "Tracking/Data/Matfiles");

% Constants
tickFontSize = 8;
labelFontSize = 10;

markerSize = 5;
turbineColor = [0.80 0.80 0.10];

% Motion-shadow appearance (oldest -> newest)
shadowSizeOldest = 1;
shadowSizeNewest = 15;
shadowAlphaOldest = 0.05;
shadowAlphaNewest = 1.00;

% Optional camera-depth scaling for the shadow markers. Scatter marker
% sizes are screen-space areas, so this explicitly adds perspective.
useShadowPerspective = true;
shadowPerspectiveStrength = 1.0; % 0 = off, 1 = physical 1/depth^2 area scaling
shadowDepthScaleLimits = [0.45, 2.25]; % [minimum, maximum] size multiplier

% Hide historical trail markers covered by the turbine's current solid
% geometry. Tests are performed in the current turbine-local coordinate
% system, so the exclusion volumes follow all six degrees of freedom. The
% newest marker at each tracked location remains visible.
useTrailExclusion = true;

% Small clearance around each analytical exclusion volume. This prevents a
% finite-size scatter marker from appearing partially embedded in a surface.
exclusionPadding = 1; % mm

% Buoys: vertical cylinders. The reference points in shadowPoints0 are at
% the centers of the buoy top surfaces.
buoyDiameter = 50; % mm
buoyDepth = 45;    % mm downward from the top surface

% Tower: elliptical cylinder from the corrected STL origin to the center of
% the hub. The major axis is assumed to be streamwise (local X), and the
% minor axis spanwise (local Y). Swap these two values if the STL uses the
% opposite orientation.
towerMajorLength = 8.4; % mm, full streamwise dimension
towerMinorLength = 6.4; % mm, full spanwise dimension
towerHeight = 100;      % mm, local Z = 0 to local Z = 100

% Hub: streamwise cylinder centered directly above the tower. Its center is
% at local [0, 0, towerHeight].
hubDiameter = 15; % mm
hubLength = 35;   % mm, cylinder axis along local X

% Number of recent frames retained in the trail
trailLength = 15;
shadowColormap = 'BuPu';

% Frames per second
frameRate = 30;   

% Seconds
duration = 5;  


% Which case to make CGI
farm_arrangement = 'Inline';
farm_spacing = 'SX50';
turbine = 2;
wave = 'LM5_AK12';
harmonic_cutoff = 2;

% Load tracking data
tracking_file = sprintf('OffshoreTracking_AllDataCombined_SavitskyGolay_Cutoff_%1.0f.mat', harmonic_cutoff);
full_tracking = load(fullfile(tracking_path, "AllData", tracking_file));
full_tracking = full_tracking.tracking;

% Load Savitsky-Golay filter design
filterspecs = full_tracking.FilterDesign;

% Load data for specific layout
tracking = full_tracking.(farm_arrangement);

% Saves
videoSavePath = '/Users/zeinsadek/Desktop/Experiments/Offshore/Tracking/Processing/TrackingCGI';
% videoSaveName = sprintf('OffshoreMotion_RealSignal_%s_%s_%s_Turbine_%s.mp4', farm_arrangement, farm_spacing, wave, num2str(turbine));
videoSaveName = sprintf('OffshoreMotion_RealSignal_%s_%s_%s_Turbine_%s_FadingGrowing_Perspective_Exclusion_Test.mp4', farm_arrangement, farm_spacing, wave, num2str(turbine));

% PNG export resolution
% Increasing this (~600) makes the videos crisper
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
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.5);

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

% Perspective projects the CAD geometry naturally. Note that scatter-marker
% sizes still require the explicit depth correction used below.
camproj(ax, 'perspective')

% Reverse MATLAB Y to match the experimental spanwise convention
set(ax, 'YDir', 'reverse')




%%% Motion-shadow points
% Each row is one reference point attached to the undeformed turbine.
% Coordinates must be in the same coordinate system as Vdata.
shadowPoints0 = [ ...
      17.5,  00.0, 100.0;   % Hub
     -41.9,  72.5, -45.0;   % Front-left buoy
     -41.9, -72.5, -45.0;   % Front-right buoy
      83.7,  00.0, -45.0;   % Rear buoy-28.868
     -25.5,  00.0, 175.0;   % Top rotor tip
     -25.5,  64.9,  62.5;   % Turbine lower left rotor tip
     -25.5, -64.9,  62.5];  % Turbine lower right rotor tip

    

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
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 'flat', ...
    'AlphaDataMapping', 'none');

colormap(ax, trailColors)




%%% Dummy six-DOF motion
nFrames   = frameRate * duration;

% Exclude duplicated endpoint
t = (0:nFrames-1) / frameRate;

%%% LOAD REAL DATA
surgeData = tracking.(farm_spacing).(wave)(turbine).x_kal * 1E3;
swayData = tracking.(farm_spacing).(wave)(turbine).y_kal * 1E3;
heaveData = tracking.(farm_spacing).(wave)(turbine).z_kal * 1E3;

rollData = tracking.(farm_spacing).(wave)(turbine).roll_kal;
pitchData = tracking.(farm_spacing).(wave)(turbine).pitch_kal;
yawData = tracking.(farm_spacing).(wave)(turbine).yaw_kal;

xlim([-100, 200])
ylim([-200, 200])
zlim([-100, 200])

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

    % Linearly vary marker size and opacity from oldest to newest.
    % During startup, use the corresponding newest portion of the full trail
    % so the current point always has the specified newest appearance.
    sizeByAge = linspace(shadowSizeOldest, shadowSizeNewest, trailLength).';
    alphaByAge = linspace(shadowAlphaOldest, shadowAlphaNewest, trailLength).';

    sizeValid = sizeByAge(end-nValid+1:end);
    alphaValid = alphaByAge(end-nValid+1:end);

    % Repeat each age value for all tracked locations
    sizePlot = repelem(sizeValid, nShadowPoints, 1);
    alphaPlot = repelem(alphaValid, nShadowPoints, 1);

    % Add camera-depth perspective to the age-based marker sizes.
    % Scatter SizeData is marker area in points^2, so true perspective-area
    % scaling varies approximately as 1/depth^2. The data-aspect-ratio
    % correction makes the camera-distance calculation robust when the axes
    % do not have equal scaling.
    if useShadowPerspective
        dar = ax.DataAspectRatio;
        Pcamera = Pplot ./ dar;
        cameraPosition = ax.CameraPosition ./ dar;
        cameraTarget = ax.CameraTarget ./ dar;

        viewDirection = cameraTarget - cameraPosition;
        viewDirection = viewDirection ./ norm(viewDirection);

        % Positive distance along the camera viewing direction
        pointDepth = (Pcamera - cameraPosition) * viewDirection.';
        referenceDepth = dot(cameraTarget - cameraPosition, viewDirection);

        % Protect against points at or behind the camera
        pointDepth = max(pointDepth, eps(referenceDepth));

        % shadowPerspectiveStrength = 0 gives no depth scaling;
        % shadowPerspectiveStrength = 1 gives 1/depth^2 area scaling.
        depthSizeScale = (referenceDepth ./ pointDepth) .^ ...
            (2 * shadowPerspectiveStrength);

        depthSizeScale = min(max(depthSizeScale, shadowDepthScaleLimits(1)), ...
                                  shadowDepthScaleLimits(2));

        sizePlot = sizePlot .* depthSizeScale;
    end

    % Hide historical trail markers that lie inside the current buoys,
    % tower, or hub. Historical world-coordinate points are transformed
    % back into the current turbine-local frame using the inverse rigid-body
    % transform. The simple analytical volumes are fast enough to evaluate
    % every animation frame and remain attached to the moving model.
    if useTrailExclusion
        Pplot_h = [Pplot, ones(size(Pplot,1),1)];
        PplotLocal_h = (M \ Pplot_h.').';
        PplotLocal = PplotLocal_h(:,1:3);

        hidePoint = false(size(PplotLocal,1),1);

        %%% Buoy exclusion: three vertical cylinders
        buoyCentersLocal = shadowPoints0(2:4,:);
        buoyRadiusExclusion = buoyDiameter/2 + exclusionPadding;

        for iBuoy = 1:size(buoyCentersLocal,1)
            center = buoyCentersLocal(iBuoy,:);

            radialDistance = hypot( ...
                PplotLocal(:,1) - center(1), ...
                PplotLocal(:,2) - center(2));

            relativeHeight = PplotLocal(:,3) - center(3);

            insideBuoy = ...
                radialDistance <= buoyRadiusExclusion & ...
                relativeHeight >= -(buoyDepth + exclusionPadding) & ...
                relativeHeight <= exclusionPadding;

            hidePoint = hidePoint | insideBuoy;
        end

        %%% Tower exclusion: elliptical vertical cylinder
        % The tower is centered on local X = 0 and Y = 0. Adding the same
        % padding to both semi-axes gives a small visual clearance around
        % the physical surface.
        towerSemiAxisX = towerMajorLength/2 + exclusionPadding;
        towerSemiAxisY = towerMinorLength/2 + exclusionPadding;

        towerEllipseCoordinate = ...
            (PplotLocal(:,1) ./ towerSemiAxisX).^2 + ...
            (PplotLocal(:,2) ./ towerSemiAxisY).^2;

        insideTower = ...
            towerEllipseCoordinate <= 1 & ...
            PplotLocal(:,3) >= -exclusionPadding & ...
            PplotLocal(:,3) <= towerHeight + exclusionPadding;

        hidePoint = hidePoint | insideTower;

        %%% Hub exclusion: finite streamwise cylinder
        % The hub cylinder is centered over the tower at local
        % [0, 0, towerHeight], with its axis along local X.
        hubCenterLocal = [0, 0, towerHeight];
        hubHalfLengthExclusion = hubLength/2 + exclusionPadding;
        hubRadiusExclusion = hubDiameter/2 + exclusionPadding;

        hubAxialDistance = abs(PplotLocal(:,1) - hubCenterLocal(1));
        hubRadialDistance = hypot( ...
            PplotLocal(:,2) - hubCenterLocal(2), ...
            PplotLocal(:,3) - hubCenterLocal(3));

        insideHub = ...
            hubAxialDistance <= hubHalfLengthExclusion & ...
            hubRadialDistance <= hubRadiusExclusion;

        hidePoint = hidePoint | insideHub;

        % Retain the newest marker for each tracked point. Only historical
        % trail markers are occluded by the current turbine geometry.
        newestPointIndices = (size(Pplot,1)-nShadowPoints+1):size(Pplot,1);
        hidePoint(newestPointIndices) = false;

        alphaPlot(hidePoint) = 0;
    end

    % Update scatter object
    hShadow.XData = Pplot(:,1);
    hShadow.YData = Pplot(:,2);
    hShadow.ZData = Pplot(:,3);
    hShadow.CData = Cplot;
    hShadow.SizeData = sizePlot;
    hShadow.AlphaData = alphaPlot;

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