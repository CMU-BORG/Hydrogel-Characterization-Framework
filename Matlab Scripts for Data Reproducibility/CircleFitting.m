%% read image

%clear all; close all; %clc

currentFolder = pwd;
root = currentFolder + "/Experimental Data/";

inventory_file = "Inventory Final.xlsx";

inventory = readtable(root+inventory_file,"ReadVariableNames",true);

A = zeros(length(inventory.Folder), 25);

%% loop through samples
for sample = 1:(length(inventory.Folder))

    if isnan(inventory.CrosslinkTime(sample))
        A(sample) = zeros(1, length(A(sample)));
        %continue
    end 

    folder = root + inventory.Folder(sample);  % access sample folder
    files = {dir(folder).name}; % read files in sample folder

    frames_folder = root + inventory.Folder(sample) + "/frames";  % access frames folder
    frame_files = {dir(frames_folder).name}; % read files in folder
    frame_files = natsort(frame_files); % fix sorting

    sample_name = inventory.SampleID(sample);        % record sample ID
    sample_name = sample_name{1}; 
    

    %% loop through frames
    t = 0;
    col = 1;
    % for every frame file in the sample folder
    for i = 1:length(frame_files)
       
        frame = frame_files{i};
        if frame((1)) ~= 't'
            continue
        end

        array = linspace(0,2*pi,100);

        img = imread(frames_folder + "/" + frame);
        hsv_img = rgb2hsv(img);

        %% find blue parts

        % Define thresholds for the blue color
        % optionally set specific thresholds for specific samples (based on lighting condition changes, reflections, etc)
        if sample == 1
            hue_min = 0.53;
            hue_max = 0.57;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3; % Adjust based on the brightness level
            
        elseif 
            hue_min = 0.556;
            hue_max = 0.58;
            sat_min = 0.95; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        elseif sample == 4 || sample == 6 
            hue_min = 0.545;
            hue_max = 0.57;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3; % Adjust based on the brightness level
        elseif sample == 5 & i < 5
            hue_min = 0.551;
            hue_max = 0.58;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        elseif sample == 5 &  i >=5 & i < 7
            hue_min = 0.5535;
            hue_max = 0.58;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        elseif sample == 5 &  i >= 7
            hue_min = 0.55341;
            hue_max = 0.58;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        elseif sample == 8 
            hue_min = 0.556;
            hue_max = 0.58;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        elseif sample == 9
            hue_min = 0.556;
            hue_max = 0.58;
            sat_min = 0.99; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        elseif sample == 11
            hue_min = 0.552;
            hue_max = 0.585;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        elseif sample == 13 & i == 8
            hue_min = 0.5415;
            hue_max = 0.57;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        elseif sample == 13
            hue_min = 0.53;
            hue_max = 0.57;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3; % Adjust based on the brightness level
        elseif sample == 21 & i <= 11
            hue_min = 0.539;
            hue_max = 0.58;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3; % Adjust based on the brightness level
        elseif sample == 30
            hue_min = 0.56;
            hue_max = 0.58;
            sat_min = 0.99; % Adjust based on the saturation level of blue in your image
            val_min = 0.5;
        else
            hue_min = 0.552;
            hue_max = 0.58;
            sat_min = 0.91; % Adjust based on the saturation level of blue in your image
            val_min = 0.3;
        end

        % Extract the individual HSV channels
        hue = hsv_img(:,:,1);
        saturation = hsv_img(:,:,2);
        value = hsv_img(:,:,3);

        
        if sample == 25 & i >=15
            img = img(1:1000, :, :);
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 10) & (R < 90) & (G < 95);  % Adjust thresholds to suit your image
        elseif sample == 25
            img = img(1:1000, :, :);
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 10) & (R < 120) & (G < 115);  % Adjust thresholds to suit your image
        elseif sample == 14
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 100) & (R < 100) & (G < 180);  % Adjust thresholds to suit your image
        elseif sample == 20
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 100) & (R < 100) & (G < 110);  % Adjust thresholds to suit your image
        elseif sample == 32
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 100) & (R < 50) & (G < 150);  % Adjust thresholds to suit your image
        elseif sample == 34
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 100) & (R < 40) & (G < 160);  % Adjust thresholds to suit your image
        elseif sample == 41 
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 150) & (R < 130) & (G < 200);  % Adjust thresholds to suit your image
        elseif sample == 44 
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 100) & (R < 90) & (G < 180);  % Adjust thresholds to suit your image
        elseif sample == 47 
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 190) & (R < 175) & (G < 200);  % Adjust thresholds to suit your image
        elseif sample >= 36  
            R = img(:,:,1);
            G = img(:,:,2);
            B = img(:,:,3);
            blue_mask = (B > 100) & (R < 100) & (G < 180);  % Adjust thresholds to suit your image
        else
            % Create a binary mask where the blue color is located
            blue_mask = (hue >= hue_min) & (hue <= hue_max) & ...
                (saturation >= sat_min) & ...
                (value >= val_min);
        end

        %% filter out small bits
        % Initialize an empty image to hold the result
        blue_isolated = zeros(size(img), 'like', img);

        % Copy the original image's blue regions into the new image
        for j = 1:3
            blue_isolated(:,:,j) = img(:,:,j) .* uint8(blue_mask);
        end

        % Define a minimum size threshold for the blue regions (in pixels)
        min_area = 10000;  % Adjust this value based on your needs

        % Remove small blue regions
        blue_mask_filtered = bwareaopen(blue_mask, min_area);

        % Optional: Smooth out the mask using morphological operations
        se = strel('disk', 5);  % Structuring element, adjust size as needed
        blue_mask_filtered = imdilate(blue_mask_filtered, se);  % Dilate to expand blue areas
        blue_mask_filtered = imerode(blue_mask_filtered, se);   % Erode to clean up small areas

        % Initialize an empty image to hold the result
        blue_isolated_filtered = zeros(size(img), 'like', img);

        % Apply the filtered mask to the original image
        for j = 1:3
            blue_isolated_filtered(:,:,j) = img(:,:,j) .* uint8(blue_mask_filtered);
        end

        %% Fit outer circles using edge detection

        diam_data = mean(blue_isolated_filtered,3); % convert to greyscale

        [x_,y_] = size(diam_data);

        thres = mean(diam_data,"all"); % threshold for black vs white based on mean
        data_BW = imgaussfilt(diam_data, 0.5) > thres; % Smooth and binarize

        % detect all edges
        edge_data = edge(data_BW);
        [xx,yy] = find(edge_data);

        % isolate inner region
        min_x = min(xx);
        max_x = max(xx);
        range_x = max_x - min_x;
        min_y = min(yy);
        max_y = max(yy);
        range_y = max_y - min_y;

        upper = floor(min_x + range_x/4);
        lower = floor(max_x - range_x/4);
        left = floor(min_y + range_y/4);
        right = floor(max_y - range_y/4);

        % create logical arrays the size of edge_data but all zeros
        inner_edge = false(size(edge_data));
        outer_edge = edge_data;

        % Set the logical array to true for the defined region of true values from
        % edge_data
        inner_edge(upper:lower, left:right) = edge_data(upper:lower, left:right);
        [xxi, yyi] = find(inner_edge);
        outer_edge(upper:lower, left:right) = 0;
        [xxo, yyo] = find(outer_edge);

        % initial guess for circle center
        x0 = mean(xx);
        y0 = mean(yy);
        [xx_bw,yy_bw] = find(data_BW);
        x0 = mean(xx_bw);
        y0 = mean(yy_bw);

        % Fit circle to the outer edge
        outer_err = @(X) ErrEllipse(xxo(:), yyo(:), X(1), X(1), X(3), X(4));
        params_outer = fminsearch(outer_err, [x_ / 8, x_ / 8, x0, y0]);
        rx_outer = params_outer(1);
        ry_outer = params_outer(1);
        x0 = params_outer(3);
        y0 = params_outer(4);
        R1 = mean([rx_outer, ry_outer]);

        % Fit circle to the inner edge
        inner_err = @(X) ErrEllipse(xxi(:), yyi(:), X(1), X(1), X(3), X(4));
        params_inner = fminsearch(inner_err, [R1/2, R1/2, x0, y0]);
        rx_inner = params_inner(1);
        ry_inner = params_inner(1);
        x0 = params_inner(3);
        y0 = params_inner(4);
        R2 = mean([rx_inner, ry_inner]);

        % fprintf("Outer Radius: %.3f px\n",R1)
        % fprintf("Inner Radius: %.3f px\n",R2)

        %% Plot the detected outer and inner circles

        process = figure; %('visible', 'off')
        subplot(2, 3, 1), imshow(img), title('Original Image');
        subplot(2, 3, 2), imshow(blue_isolated), title('Isolated Blue (Before Filtering)');
        subplot(2, 3, 3), imshow(blue_isolated_filtered), title('Isolated Blue (After Filtering)');
        subplot(2, 3, 4), imshow(diam_data), title('BW Image');
        subplot(2, 3, 5), imshow(edge_data), title('Edge Detection');
        subplot(2, 3, 6); hold on
        imshow((diam_data),'XData',[0,y_],'YData',[0,x_]);
        set(gca, 'YDir', 'normal');
        plot(yy,xx,'.g',"MarkerSize",2);
        title('Edges Overlayed');

        %% for each frame, record surface area and traced circles

        % some samples have scale shifts from video bumps
        if sample == 2
            if col > 1
                scale = inventory.Scale_pxpmm_(sample) / 0.99;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 5
            if col >= 5
                scale = 20.9;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 9
            if col > 1
                scale = inventory.Scale_pxpmm_(sample) / 1.02;
            else 
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 11
            if col >= 10
                scale = inventory.Scale_pxpmm_(sample) / 0.99;
            elseif col > 3
                scale = inventory.Scale_pxpmm_(sample) / 1.04;
            else 
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 12
            if col >= 2
                scale = inventory.Scale_pxpmm_(sample) / 0.96;
            else 
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 13
            if col >= 9
                scale = inventory.Scale_pxpmm_(sample) / 1.035;
            elseif col >= 7
                scale = inventory.Scale_pxpmm_(sample) / 1.023;
            else 
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 15
            if col >= 9
                scale = 13.223;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end 
        elseif sample == 17
            if col >= 10
                scale = inventory.Scale_pxpmm_(sample) / 0.986;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end 
        elseif sample == 19
            if col >= 10
                scale = inventory.Scale_pxpmm_(sample) / 0.988;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end 
        elseif sample == 23
            if col < 13
                scale = inventory.Scale_pxpmm_(sample);
            elseif col >= 13
                scale = inventory.Scale_pxpmm_(sample) / 0.995;
            end
        elseif sample == 24
            if col >= 13
                scale = inventory.Scale_pxpmm_(sample) / 0.985;
            elseif col >= 6
                scale = inventory.Scale_pxpmm_(sample);
            elseif col >= 5
                scale = inventory.Scale_pxpmm_(sample) / 1.01;
            elseif col >= 4
                scale = inventory.Scale_pxpmm_(sample) / 1.015;
            elseif col >= 3
                scale = inventory.Scale_pxpmm_(sample) / 1.02;
            elseif col >= 2
                scale = inventory.Scale_pxpmm_(sample) / 1.03;
            elseif col >= 1
                scale = inventory.Scale_pxpmm_(sample) / 1.06;
            end
        elseif sample == 25
            if col >= 13
                scale = inventory.Scale_pxpmm_(sample) / 1.3;
            elseif col >= 9
                scale = inventory.Scale_pxpmm_(sample) / 0.99;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 27
            if col >= 2
                scale = inventory.Scale_pxpmm_(sample);
            else 
                scale = 14.677;
            end 
        elseif sample == 28
            if col < 9 
                scale = inventory.Scale_pxpmm_(sample);
            else 
                scale = inventory.Scale_pxpmm_(sample) / 0.915;
            end
        elseif sample == 30
            if col < 9
                scale = inventory.Scale_pxpmm_(sample);
            else
                scale = inventory.Scale_pxpmm_(sample) / 0.96;
            end
        elseif sample == 40
            if col < 9
                scale = inventory.Scale_pxpmm_(sample);
            else
                scale = inventory.Scale_pxpmm_(sample) / 1.014;
            end
        elseif sample == 41
            if col >= 14
                scale = inventory.Scale_pxpmm_(sample) / 1.055;
            elseif col >=9
                scale = inventory.Scale_pxpmm_(sample) / 1.045;
            elseif col >= 8
                scale = inventory.Scale_pxpmm_(sample) / 1.005;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 42
            if col < 2
                scale = 21.186;
            elseif col >=9
                scale = inventory.Scale_pxpmm_(sample) / 0.94;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 43
            if col >=18
                scale = inventory.Scale_pxpmm_(sample) / .993;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 44
            if col >= 14
                scale = inventory.Scale_pxpmm_(sample) / 0.985;
            elseif col >= 9
                scale = inventory.Scale_pxpmm_(sample) / 1.01;
            elseif col >= 2
                scale = inventory.Scale_pxpmm_(sample) / 0.97;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 45
            if col >=22 
                scale = inventory.Scale_pxpmm_(sample) / 0.945;
            elseif col >= 18
                scale = inventory.Scale_pxpmm_(sample) / 0.95;
            elseif col >= 14
                scale = inventory.Scale_pxpmm_(sample) / 0.99;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 46
            if col >=3
                scale = inventory.Scale_pxpmm_(sample) / 0.96;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        elseif sample == 47
            if col >=2
                scale = inventory.Scale_pxpmm_(sample) / 0.94;
            else
                scale = inventory.Scale_pxpmm_(sample);
            end
        else
            scale = inventory.Scale_pxpmm_(sample);
        end

        % convert radii from px to mm
        R1_mm = R1 / scale;
        R2_mm = R2 / scale;

        if col == 1
            R1i = R1;
            R2i = R2;
        end

        row = sample;
        row;
        col;
  
        % calculate area
        A(row, col) = (pi * R1_mm^2) - (pi * R2_mm^2);

        % plot circles
        trace = figure('visible', 'off'); %('visible', 'off')
        hold on

        plot(xx,yy,'.k');
        plot(x0 + R1*cos(array), y0 + R1*sin(array),'-r',"LineWidth",2);
        plot(x0 + R2*cos(array), y0 + R2*sin(array),'-r',"LineWidth",2);
        plot(x0 + R1i*cos(array), y0 + R1i*sin(array),'--r',"LineWidth",1);
        plot(x0 + R2i*cos(array), y0 + R2i*sin(array),'--r',"LineWidth",1);
        title('Sample Radius Detection');
        ylim([y0 - 175, y0 + 175]);
        %xlim([x0 - 200, x0 + 200]);
        xticks([-10, -5, 0, 5, 10]);
        yticks([-10, -5, 0, 5, 10]);
        xlabel("px");
        ylabel("px");
        str = {"Area: " + round(A(row, col), 2) + "mm^2"};
        text(x0 - 190, y0 + 150, str);

        axis equal
        %saveas(trace,folder + "/figs/t" + t + ".png",'png'); % save svgs to frame folder

        t = t + 15;
        col = col + 1;

   end % end to looping through frames

   sample
   % save filtering and masking figure of the last frame
   %saveas(process,folder + "/figs/Process.svg",'svg'); % save svgs to drive folder for last frame

   % %% HSV Mapping of last frame
   % % Extract the V (value) channel from the HSV array
   % V_channel = blue_isolated_filtered(:,:,3);
   % 
   % % Extract the intensity values along the row y0
   % intensity_values = V_channel(floor(x0), (min(yyo)-50):(max(yyo)+50));
   % 
   % % Create the plot
   % hsvMap = figure('visible', 'off');
   % plot(intensity_values);
   % xlabel('Column Index');
   % ylabel('Intensity');
   % title(sprintf('Intensity of HSV Pixels Along Row %d', x0));
   % grid on;
   % saveas(hsvMap,folder + "/figs/hsvMap.svg",'svg'); % save hsv svgs to drive folder for last frame

 

%end % end to looping through samples

%writematrix(A, 'Areas Final new.csv');


