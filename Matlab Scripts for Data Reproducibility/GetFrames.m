clear all; close all; %clc

currentFolder = pwd;
root = currentFolder + "/Experimental Data/";

inventory_file = "Inventory.xlsx";

inventory = readtable(root+inventory_file,"ReadVariableNames",true);

%% read inventory and select video file

sample = 25;
    folder = root + inventory.Folder(sample) + "/";  % access sample folder
    files = {dir(folder).name};                 % read files in folder
    videofile = files{4};                         % select video
    
    sample_name = inventory.SampleID(sample);        % record sample ID
    sample_name = sample_name{1};               % not sure what this is doing

    video = VideoReader(folder + videofile);
    numFrames = video.FrameRate * video.Duration;
    interval = video.Framerate * 900; % frames per 15 min intervals
    %numintervals = floor(numFrames / interval); % number of 15 min segments

    %% get frames every 15 min
    j = 1;
    t = 0;

    while j <= numFrames
        frame = read(video, int64(j));
        time = "t" + t + ".jpg";
        imwrite(frame, folder + "frames/" + time);
        j = j + interval
        t = t + 15
    end

