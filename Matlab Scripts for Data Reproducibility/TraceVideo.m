%% trace video 

clear all; close all; %clc

% Define folder path containing images (converted from SVGs)
root = "/Users/averywilliamson/Library/CloudStorage/GoogleDrive-averyw@andrew.cmu.edu/Shared drives/CMU BORG - Hydrogel Printing/ALGINATE CHARACTERIZATION/MTS Experiments/Ring/";

inventory_file = "Inventory.xlsx";

inventory = readtable(root+inventory_file,"ReadVariableNames",true);

imageFolder = root + inventory.Folder(21) + "/figs";

imageFiles = dir(fullfile(imageFolder, '*.png')); % Replace with *.jpg if you have JPEGs

% Sort files if necessary
imageFiles = natsort({imageFiles.name}); % Optional: if you want to sort the images numerically

% Set up the video writer
outputVideo = VideoWriter('fade_through_svgs', 'MPEG-4'); % Specify the output file name
outputVideo.FrameRate = 10; % Adjust the frame rate as needed
open(outputVideo);

% speed = 9000x (2.5 hours per second)

% Read and store the images
numImages = length(imageFiles);
imageSequence = cell(numImages, 1);

for i = 1:numImages
    img = imread(fullfile(imageFolder, imageFiles{i}));
    imageSequence{i} = im2double(img); % Convert to double for smooth blending
end

% Define the fade steps between images (e.g., 10 steps per transition)
fadeSteps = 10;

% Loop through each image and create fade transitions
for i = 1:numImages-1
    img1 = imageSequence{i};
    img2 = imageSequence{i+1};
    
    % Create fade effect between img1 and img2
    for alpha = linspace(0, 1, fadeSteps)
        blendedImg = (1 - alpha) * img1 + alpha * img2;
        writeVideo(outputVideo, blendedImg);
    end
end

% Write the final image to the video to hold the last frame for a moment
for holdFrame = 1:fadeSteps
    writeVideo(outputVideo, imageSequence{end});
end

% Close the video file
close(outputVideo);

disp('Video has been created successfully!');