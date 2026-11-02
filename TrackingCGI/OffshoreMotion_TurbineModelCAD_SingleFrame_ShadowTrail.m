%% Single-frame turbine motion rendering with a historical shadow trail
% Renders one selected measured-motion frame together with a fading and
% shrinking trail from a user-specified number of preceding frames.

clc; close all; clear

% Paths
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
offshore_path = "/Users/zeinsadek/Desktop/Experiments/Offshore";
tracking_path = fullfile(offshore_path, "Tracking/Data/Matfiles");

% Constants
tickFontSize = 8;
labelFontSize = 10;

markerSize = 5;
% turbineColor = [0.80 0.80 0.10];
turbineColor = [0.42, 0.44, 0.46];

% Motion-shadow appearance (oldest -> newest)
shadowSizeOldest = 1;
shadowSizeNewest = 20;
shadowAlphaOldest = 0.2;
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
useTrailExclusion = false;

% Small clearance around each analytical exclusion volume. This prevents a
% finite-size scatter marker from appearing partially embedded in a surface.
exclusionPadding = 2; % mm

% Buoys: vertical cylinders. Their top-center locations are defined
% independently below in buoyCenters0.
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

% Independent exclusion-volume locations in the undeformed turbine-local
% coordinate system. These are intentionally separate from shadowPoints0,
% so shadow points can be added, removed, or reordered freely.
%
% The three buoy centers form an equilateral triangle with 100 mm side
% length. Coordinates are the centers of the buoy top surfaces.
buoyCenters0 = [ ...
    -28.868,  50.000, 0.0; ...
    -28.868, -50.000, 0.0; ...
     57.735,   0.000, 0.0];

% Tower cross-section center and hub-cylinder center.
towerCenter0 = [0, 0];
hubCenter0 = [0, 0, towerHeight];

% Single-image controls
% targetFrame is the measured-data sample shown by the solid CAD model.
% priorFrames is the number of samples before targetFrame included in the
% trail. The target frame itself is also plotted, so the total number of
% trail ages is priorFrames + 1.
targetFrame = 360;
priorFrames = 15;
% shadowColormap = 'amethyst';
shadowColormap = 'wistia';

% Sampling rate, used only to report the selected time in the title/name.
frameRate = 30;

% 3D viewing angle
viewAz = -21;
viewEl = 40;




%%% Which case to make CGI
farm_arrangement = 'Staggered';
farm_spacing = 'SX50';
turbine = 5;
wave = 'LM5_AK12';
harmonic_cutoff = 5;

% Load tracking data
tracking_file = sprintf('OffshoreTracking_AllDataCombined_SavitskyGolay_Cutoff_%1.0f.mat', harmonic_cutoff);
full_tracking = load(fullfile(tracking_path, "AllData", tracking_file));
full_tracking = full_tracking.tracking;

% Load Savitsky-Golay filter design
filterspecs = full_tracking.FilterDesign;

% Load data for specific layout
tracking = full_tracking.(farm_arrangement);

% Save location
caseFolder = [wave, '_', farm_spacing, '_', farm_arrangement];
imageSavePath = fullfile( ...
    '/Users/zeinsadek/Desktop/Experiments/Offshore/Tracking/Processing/TrackingCGI/SingleImages', ...
    caseFolder);
if ~exist(imageSavePath, 'dir'); mkdir(imageSavePath); end

% PNG export resolution
frameResolution = 600;


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
view(ax, viewAz, viewEl)
grid(ax, 'on')

camlight(ax, 'headlight')
material(ax, 'dull')

% Perspective projects the CAD geometry naturally. Note that scatter-marker
% sizes still require the explicit depth correction used below.
camproj(ax, 'perspective')

% Reverse MATLAB Y to match the experimental spanwise convention
set(ax, 'YDir', 'reverse')




%%% Motion-shadow points
% Each row is one reference point attached to the undeformed turbine.
% Coordinates must be in the same coordinate system as Vdata. These points
% are independent of the exclusion geometry, so they can be freely added,
% removed, or reordered without changing the buoy/tower/hub masks.
% shadowPoints0 = [ ...
%       17.5,  00.0, 100.0;   % Hub
%      -41.9,  72.5,   0.0;   % Front-left buoy
%      -41.9, -72.5,   -45;   % Front-right buoy
%       83.7,  00.0,   -45;   % Rear buoy
%      -25.5,  00.0, 175.0;   % Top rotor tip
%      -25.5,  64.9,  62.5;   % Turbine lower left rotor tip
%      -25.5, -64.9,  62.5];  % Turbine lower right rotor tip

% shadowPoints0 = [ ...
%       % 17.5,  00.0, 100.0;   % Hub
%      -41.9,  72.5,   0.0;   % Front-left buoy
%      -41.9, -72.5,   0.0;   % Front-right buoy
%       83.7,  00.0,   0.0];   % Rear buoy
%      % -25.5,  00.0, 175.0;   % Top rotor tip
%      % -25.5,  64.9,  62.5;   % Turbine lower left rotor tip
%      % -25.5, -64.9,  62.5];  % Turbine lower right rotor tip

%%% Parameterized side length
height = 0;
sideLength = 153;   % Triangle side length

% Height of the equilateral triangle
triangleHeight = sqrt(3)/2 * sideLength;

% Corner coordinates [x, y]
shadowPoints0 = [ ...
     2*triangleHeight/3, 0, height;   % Top corner
    -triangleHeight/3, -sideLength/2, height;     % Bottom-left corner
    -triangleHeight/3,  sideLength/2,  height];    % Bottom-right corner
    

nShadowPoints = size(shadowPoints0, 1);

% Total number of retained ages, including the selected current frame
totalTrailFrames = priorFrames + 1;

% Dimensions: totalTrailFrames x nShadowPoints x 3
pointHistory = nan(totalTrailFrames, nShadowPoints, 3);

% Colormap for oldest-to-newest points
trailColors = flipud(slanCM(shadowColormap, totalTrailFrames));

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




%%% LOAD REAL DATA
surgeData = tracking.(farm_spacing).(wave)(turbine).x_kal * 1E3;
swayData = tracking.(farm_spacing).(wave)(turbine).y_kal * 1E3;
heaveData = tracking.(farm_spacing).(wave)(turbine).z_kal * 1E3;

rollData = tracking.(farm_spacing).(wave)(turbine).roll_kal;
pitchData = tracking.(farm_spacing).(wave)(turbine).pitch_kal;
yawData = tracking.(farm_spacing).(wave)(turbine).yaw_kal;







% Plot limits
xlim(ax, [-100, 150])
ylim(ax, [-200, 100])
zlim(ax, [-100, 175])
axis(ax, 'manual')

% Axes labels
% xlabel(ax, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel(ax, '$z$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% zlabel(ax, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Axes tick markers
xticks(-50:50:150)
yticks(-200:100:100)
zticks(-50:50:150)


%% Validate the requested frame and assemble its history
nAvailable = min([ ...
    numel(surgeData), numel(swayData), numel(heaveData), ...
    numel(rollData), numel(pitchData), numel(yawData)]);

if targetFrame < 1 || targetFrame > nAvailable || targetFrame ~= round(targetFrame)
    error('targetFrame must be an integer from 1 to %d.', nAvailable);
end

firstTrailFrame = max(1, targetFrame - priorFrames);
trailFrameIndices = firstTrailFrame:targetFrame;
nValid = numel(trailFrameIndices);

if nValid < totalTrailFrames
    warning(['Only %d prior frames are available before targetFrame = %d. ' ...
        'The trail will contain %d total frames.'], ...
        targetFrame - 1, targetFrame, nValid);
end

% Use the matching newest portion of each appearance ramp when the target
% frame is close to the start of the measured record.
Cvalid = trailColors(end-nValid+1:end, :);
sizeByAge = linspace(shadowSizeOldest, shadowSizeNewest, totalTrailFrames).';
alphaByAge = linspace(shadowAlphaOldest, shadowAlphaNewest, totalTrailFrames).';
sizeValid = sizeByAge(end-nValid+1:end);
alphaValid = alphaByAge(end-nValid+1:end);

%% Transform the tracked reference points at each retained frame
rotationOrigin = [0, 0, 0];

T_toOrigin = eye(4);
T_toOrigin(1:3,4) = -rotationOrigin(:);

T_fromOrigin = eye(4);
T_fromOrigin(1:3,4) = rotationOrigin(:);

P0_h = [shadowPoints0, ones(nShadowPoints,1)];
Mcurrent = eye(4);

for k = 1:nValid
    n = trailFrameIndices(k);

    % Translations [mm]
    surge = surgeData(n);       % MATLAB X = experimental x
    sway  = -swayData(n);       % MATLAB Y = experimental z
    heave = heaveData(n);       % MATLAB Z = experimental y

    % Rotations [rad]
    roll  = -deg2rad(rollData(n));
    pitch = -deg2rad(pitchData(n));
    yaw   = -deg2rad(yawData(n));

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

    T_motion = eye(4);
    T_motion(1:3,4) = [surge; sway; heave];

    M = T_motion * T_fromOrigin * R * T_toOrigin;

    P_h = (M * P0_h.').';
    pointHistory(k,:,:) = reshape(P_h(:,1:3), 1, nShadowPoints, 3);

    if n == targetFrame
        Mcurrent = M;
    end
end

% Place the solid CAD model at the selected current frame.
motionTransform.Matrix = Mcurrent;

%% Build the complete shadow scatter data
% Flatten as age 1/all points, age 2/all points, etc.
historyValid = pointHistory(1:nValid,:,:);
Pplot = reshape(permute(historyValid, [2, 1, 3]), [], 3);

Cplot = repelem(Cvalid, nShadowPoints, 1);
sizePlot = repelem(sizeValid, nShadowPoints, 1);
alphaPlot = repelem(alphaValid, nShadowPoints, 1);

% Add camera-depth perspective to the age-based marker sizes.
if useShadowPerspective
    dar = ax.DataAspectRatio;
    Pcamera = Pplot ./ dar;
    cameraPosition = ax.CameraPosition ./ dar;
    cameraTarget = ax.CameraTarget ./ dar;

    viewDirection = cameraTarget - cameraPosition;
    viewDirection = viewDirection ./ norm(viewDirection);

    pointDepth = (Pcamera - cameraPosition) * viewDirection.';
    referenceDepth = dot(cameraTarget - cameraPosition, viewDirection);
    pointDepth = max(pointDepth, eps(referenceDepth));

    depthSizeScale = (referenceDepth ./ pointDepth) .^ ...
        (2 * shadowPerspectiveStrength);

    depthSizeScale = min(max(depthSizeScale, shadowDepthScaleLimits(1)), ...
                              shadowDepthScaleLimits(2));

    sizePlot = sizePlot .* depthSizeScale;
end

%% Hide historical markers covered by the current solid geometry
if useTrailExclusion
    Pplot_h = [Pplot, ones(size(Pplot,1),1)];
    PplotLocal_h = (Mcurrent \ Pplot_h.').';
    PplotLocal = PplotLocal_h(:,1:3);

    hidePoint = false(size(PplotLocal,1),1);

    % Buoys: vertical cylinders
    buoyRadiusExclusion = buoyDiameter/2 + exclusionPadding;
    for iBuoy = 1:size(buoyCenters0,1)
        center = buoyCenters0(iBuoy,:);

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

    % Tower: elliptical vertical cylinder
    towerSemiAxisX = towerMajorLength/2 + exclusionPadding;
    towerSemiAxisY = towerMinorLength/2 + exclusionPadding;

    towerEllipseCoordinate = ...
        ((PplotLocal(:,1) - towerCenter0(1)) ./ towerSemiAxisX).^2 + ...
        ((PplotLocal(:,2) - towerCenter0(2)) ./ towerSemiAxisY).^2;

    insideTower = ...
        towerEllipseCoordinate <= 1 & ...
        PplotLocal(:,3) >= -exclusionPadding & ...
        PplotLocal(:,3) <= towerHeight + exclusionPadding;

    hidePoint = hidePoint | insideTower;

    % Hub: finite streamwise cylinder
    hubHalfLengthExclusion = hubLength/2 + exclusionPadding;
    hubRadiusExclusion = hubDiameter/2 + exclusionPadding;

    hubAxialDistance = abs(PplotLocal(:,1) - hubCenter0(1));
    hubRadialDistance = hypot( ...
        PplotLocal(:,2) - hubCenter0(2), ...
        PplotLocal(:,3) - hubCenter0(3));

    insideHub = ...
        hubAxialDistance <= hubHalfLengthExclusion & ...
        hubRadialDistance <= hubRadiusExclusion;

    hidePoint = hidePoint | insideHub;

    % Always retain the newest/current marker at every tracked location.
    newestPointIndices = (size(Pplot,1)-nShadowPoints+1):size(Pplot,1);
    hidePoint(newestPointIndices) = false;
    alphaPlot(hidePoint) = 0;
end

%% Draw and export the selected image
hShadow.XData = Pplot(:,1);
hShadow.YData = Pplot(:,2);
hShadow.ZData = Pplot(:,3);
hShadow.CData = Cplot;
hShadow.SizeData = sizePlot;
hShadow.AlphaData = alphaPlot;

selectedTime = (targetFrame - 1) / frameRate;
% title(ax, sprintf('$t = %.2f$ s', selectedTime), 'Interpreter', 'latex');

drawnow

imageSaveName = sprintf( ...
    'OffshoreMotion_%s_%s_%s_Turbine_%d_Frame_%d_PriorFrames_%d.png', ...
    farm_arrangement, farm_spacing, wave, turbine, targetFrame, priorFrames);
imageFile = fullfile(imageSavePath, imageSaveName);

exportgraphics(fig, imageFile, ...
    'Resolution', frameResolution, ...
    'BackgroundColor', 'white');

fprintf('Single-frame rendering saved successfully:\n%s\n', imageFile);
fprintf('Current frame: %d (t = %.3f s)\n', targetFrame, selectedTime);
fprintf('Trail frames used: %d through %d (%d total)\n', ...
    firstTrailFrame, targetFrame, nValid);
