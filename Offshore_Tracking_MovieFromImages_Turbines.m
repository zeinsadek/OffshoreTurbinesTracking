%% Create video from TIFF image sequence
% Converts all .tif/.tiff images in a folder into an MPEG-4 video.

clear;
clc;

%% User inputs

arrangement = 'Staggered';
% HDD = 'E:Tracking';
HDD = '/Volumes/FWF_SX50_S/Tracking';
HDD = [HDD, '_', arrangement];
spacing = 'WT60_SX50_AG0';
wave = 'LM5_AK12';
camera = 'CAM4';

harmonic_cases = {'LM25_AK06', 'LM25_AK09', 'LM25_AK12', ...
                  'LM33_AK06', 'LM33_AK09', 'LM33_AK12'};

if strcmp(arrangement, 'Inline')
    if any(strcmp(wave, harmonic_cases))
        case_folder = [spacing, '_Harmonics'];
    else
        case_folder = spacing;
    end
else
    case_folder = spacing;
end

inputFolder = fullfile(HDD, case_folder, wave, camera);


% Output folder/name
% outputFolder = 'C:\Users\sadek\Desktop\ZeinScratch\TrackingMovies\TurbineMovies';
outputFolder = '/Users/zeinsadek/Downloads';
outputFileName = strcat(spacing, '_', arrangement, '_', wave, '_', camera, '_TrackingMovie.mp4');

frameRate = 30;       % Video frame rate [frames/second]
videoQuality = 100;   % MPEG-4 quality from 0 to 100

% Crop frame dimensions to a multiple of this value.
% A multiple of 8 is generally safe for MPEG-4 encoding.
dimensionMultiple = 8;

%% Find TIFF images

filesTif  = dir(fullfile(inputFolder, '*.tif'));
filesTiff = dir(fullfile(inputFolder, '*.tiff'));

imageFiles = [filesTif; filesTiff];

if isempty(imageFiles)
    error('No .tif or .tiff images were found in:\n%s', inputFolder);
end

%% Sort images by filename

fileNames = {imageFiles.name};

% Natural sorting when natsortfiles is available.
% For example: frame2.tif comes before frame10.tif.
if exist('natsortfiles', 'file') == 2
    [~, sortIndex] = natsortfiles(fileNames);
else
    % MATLAB's normal alphabetical sorting
    [~, sortIndex] = sort(lower(fileNames));

    warning(['natsortfiles was not found. Files will be sorted ', ...
             'alphabetically rather than numerically.']);
end

imageFiles = imageFiles(sortIndex);

%% Create output folder if needed

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

outputPath = fullfile(outputFolder, outputFileName);

%% Read first frame and determine video dimensions

firstImagePath = fullfile(inputFolder, imageFiles(1).name);
firstImage = imread(firstImagePath);

[firstImage, targetHeight, targetWidth] = prepareFrame( ...
    firstImage, dimensionMultiple);

fprintf('Number of frames: %d\n', numel(imageFiles));
fprintf('Original first-frame size: %d x %d pixels\n', ...
    size(imread(firstImagePath), 2), size(imread(firstImagePath), 1));
fprintf('Video frame size: %d x %d pixels\n', ...
    targetWidth, targetHeight);
fprintf('Frame rate: %.2f fps\n', frameRate);

%% Initialize VideoWriter

videoWriter = VideoWriter(outputPath, 'MPEG-4');
videoWriter.FrameRate = frameRate;
videoWriter.Quality = videoQuality;

open(videoWriter);

%% Add images to video

try
    for k = 1:numel(imageFiles)

        imagePath = fullfile(inputFolder, imageFiles(k).name);
        img = imread(imagePath);

        % Convert indexed images to RGB when necessary
        if ndims(img) == 2
            img = repmat(img, 1, 1, 3);
        elseif size(img, 3) == 4
            % Discard alpha channel
            img = img(:, :, 1:3);
        end

        % Confirm that the current image is large enough
        if size(img, 1) < targetHeight || size(img, 2) < targetWidth
            error(['Image "%s" is smaller than the required video size.\n', ...
                   'Image size: %d x %d\nRequired size: %d x %d'], ...
                   imageFiles(k).name, ...
                   size(img, 2), size(img, 1), ...
                   targetWidth, targetHeight);
        end

        % Crop every frame to exactly the same dimensions
        img = img(1:targetHeight, 1:targetWidth, :);

        writeVideo(videoWriter, img);

        fprintf('Added frame %d of %d: %s\n', ...
            k, numel(imageFiles), imageFiles(k).name);
    end

    close(videoWriter);

catch ME
    % Ensure the video file is closed if an error occurs
    close(videoWriter);
    rethrow(ME);
end

fprintf('\nVideo successfully saved to:\n%s\n', outputPath);

%% Local function

function [img, targetHeight, targetWidth] = prepareFrame(img, multiple)
%PREPAREFRAME Convert an image to RGB and crop its dimensions so that they
%are integer multiples of the specified value.

    % Convert grayscale to RGB
    if ndims(img) == 2
        img = repmat(img, 1, 1, 3);
    elseif size(img, 3) == 4
        % Remove alpha channel
        img = img(:, :, 1:3);
    end

    imageHeight = size(img, 1);
    imageWidth  = size(img, 2);

    targetHeight = imageHeight - mod(imageHeight, multiple);
    targetWidth  = imageWidth  - mod(imageWidth, multiple);

    if targetHeight <= 0 || targetWidth <= 0
        error('The image dimensions are too small for the selected multiple.');
    end

    img = img(1:targetHeight, 1:targetWidth, :);
end