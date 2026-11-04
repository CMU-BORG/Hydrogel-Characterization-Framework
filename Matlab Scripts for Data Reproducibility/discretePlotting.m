%% LOAD DATA 
clear all; close all; %clc

currentFolder = pwd;
root = currentFolder + "/Experimental Data/";

inventory_file = "Inventory Final.xlsx";
inventory = readtable(root + inventory_file,"ReadVariableNames",true);


% Extract the relevant columns from the table
crosslinkTemp = inventory.CrosslinkTemp;       % Should contain '27°C' or '37°C'
crosslinkTime = inventory.CrosslinkTime;       % Should contain numeric values like 2, 4, or 6

% Convert the CrosslinkTemp to a categorical variable for grouping by color, remove skipped samples and categorize data
crosslinkTemp = categorical(crosslinkTemp);
crosslinkTime = categorical(crosslinkTime);

unusedTime = isundefined(crosslinkTime);
inventory(unusedTime, :) = [];
crosslinkTime(unusedTime) = [];
crosslinkTime = removecats(crosslinkTime);

unusedTemp = isundefined(crosslinkTemp);
crosslinkTemp(unusedTemp) = [];
crosslinkTemp = removecats(crosslinkTemp);

uniqueTimes = unique(crosslinkTime);
uniqueTemps = unique(crosslinkTemp);

% crosslinkTemp2 = categorical(crosslinkTemp2);
% crosslinkTime2 = categorical(crosslinkTime2);
% uniqueTimes2 = unique(crosslinkTime2);
% uniqueTemps2 = unique(crosslinkTemp2);

modulusData = inventory.EffectiveElasticModulus * 1000;  % Data to be plotted
% modulusData2 = inventory2.EffectiveElasticModulus * 1000;  % Data to be plotted

initialOD = inventory.initialOD_px_ ./ inventory.Scale_pxpmm_;
initialID = inventory.initialID_px_ ./ inventory.Scale_pxpmm_;

finalOD = inventory.finalOD_px_ ./ inventory.Scale_pxpmm_;
finalID = inventory.finalID_px_ ./ inventory.Scale_pxpmm_;

initialA = (pi()*(initialOD/2).^2) - (pi()*(initialID/2).^2);
finalA = (pi()*(finalOD/2).^2) - (pi()*(finalID/2).^2);
meanInitialA = groupsummary(initialA, [crosslinkTime, crosslinkTemp],'mean');
stdInitialA = groupsummary(initialA, [crosslinkTime, crosslinkTemp],'std');
meanA0 = [meanInitialA(1)];
meanA27 = [meanInitialA(2); meanInitialA(4); meanInitialA(6)];
meanA37 = [meanInitialA(3); meanInitialA(5); meanInitialA(7)];

initialM = inventory.initialMass_g_;
finalM = inventory.finalMass_g_;
meanInitialM = groupsummary(initialM, [crosslinkTime, crosslinkTemp],'mean');
stdInitialM = groupsummary(initialM, [crosslinkTime, crosslinkTemp],'std');
meanM0 = [meanInitialM(1)];
meanM27 = [meanInitialM(2); meanInitialM(4); meanInitialM(6)];
meanM37 = [meanInitialM(3); meanInitialM(5); meanInitialM(7)];

meanModulus = groupsummary(modulusData, [crosslinkTime, crosslinkTemp],'mean');
meanE0 = [meanModulus(1)];
meanE27 = [meanModulus(2); meanModulus(4); meanModulus(6)];
meanE37 = [meanModulus(3); meanModulus(5); meanModulus(7)];

percentM = (initialM - finalM) ./ initialM * 100;
meanPercentM = groupsummary(percentM, [crosslinkTime, crosslinkTemp],'mean');
meanPM0 = [meanPercentM(1)];
meanPM27 = [meanPercentM(2); meanPercentM(4); meanPercentM(6)];
meanPM37 = [meanPercentM(3); meanPercentM(5); meanPercentM(7)];

%% INITIALIZE PLOTTING ARRAYS 

% Convert crosslinkTime to numeric values for plotting
timeNumeric = zeros(size(crosslinkTime));

% Define positions for the box plots as before
positions = [1, 2, 3, 4]; % Corresponding to '2h', '4h', and '6h'
temperatureOffset = [0, -0.15, 0.15]; % Offset for '27°C' and '37°C' to avoid overlap
temperatureOffsetBox = [0.055, -0.15, 0.0972];
colors = validatecolor(["#000000", "#4DBEEE", "#D95319"], 'multiple');

% Map the crosslink times to their respective positions
for i = 1:length(uniqueTimes)
    timeNumeric(crosslinkTime == uniqueTimes(i)) = positions(i);
end

% Adjust x-positions for each temperature condition
xShifted = timeNumeric;
xShiftedBox = timeNumeric;

for j = 1:length(temperatureOffset)
    tempCat = categories(crosslinkTemp);
    tempCat = tempCat{j}; % '27°C' or '37°C'
    xShifted(crosslinkTemp == tempCat) = xShifted(crosslinkTemp == tempCat) + temperatureOffset(j);
    xShiftedBox(crosslinkTemp == tempCat) = xShiftedBox(crosslinkTemp == tempCat) + temperatureOffsetBox(j);
end


%% STATISTICAL ANALYSIS - STIFFNESS 

fprintf("\n=== STATISTICAL ANALYSIS (Stiffness) ===\n");

numTimes = numel(uniqueTimes);
numTemps = numel(uniqueTemps);

% Containers for normality results
normalityMatrix_stiffness = false(numTimes, numTemps);

% Loop through groups using existing categories
for i = 1:numTimes
    for j = 1:numTemps

        % Logical index for this group
        idx = (crosslinkTime == uniqueTimes(i)) & ...
              (crosslinkTemp == uniqueTemps(j));

        dataTemp = modulusData(idx);

        if isempty(dataTemp)
            continue
        end

        % Shapiro–Wilk normality test
        [pSW, ~] = swft(dataTemp);
        normalityMatrix_stiffness(i,j) = pSW{2,7} > 0.05;

    end
end

% Decide global test
allNormal = all(normalityMatrix_stiffness(:));
if allNormal
    fprintf("\nNormal groups detected .\n");
else
    fprintf("\nNon‑normal groups detected > Running KRUSKAL–WALLIS.\n");
end

% PAIRWISE KRUSKAL–WALLIS BETWEEN ALL GROUPS
fprintf("\n=== Pairwise Kruskal–Wallis Tests (Alginate Crosslinking) ===\n");

% ============================
%  PAIRWISE KRUSKAL–WALLIS TABLE OUTPUT
%  ============================

% 0 hr (uncrosslinked) — no temperature distinction
E_0h = modulusData(inventory.CrosslinkTime == 0);

% 2 hr groups
E_2h_27 = modulusData(inventory.CrosslinkTime == 2 & inventory.CrosslinkTemp == 27);
E_2h_37 = modulusData(inventory.CrosslinkTime == 2 & inventory.CrosslinkTemp == 37);

% 4 hr groups
E_4h_27 = modulusData(inventory.CrosslinkTime == 4 & inventory.CrosslinkTemp == 27);
E_4h_37 = modulusData(inventory.CrosslinkTime == 4 & inventory.CrosslinkTemp == 37);

% 6 hr groups
E_6h_27 = modulusData(inventory.CrosslinkTime == 6 & inventory.CrosslinkTemp == 27);
E_6h_37 = modulusData(inventory.CrosslinkTime == 6 & inventory.CrosslinkTemp == 37);   % <-- now included

% REGISTER GROUPS
groupNames = [
    "0h"
    "2h_27"
    "2h_37"
    "4h_27"
    "4h_37"
    "6h_27"
    "6h_37"
];

groupData = {
    E_0h
    E_2h_27
    E_2h_37
    E_4h_27
    E_4h_37
    E_6h_27
    E_6h_37
};

% PAIRWISE KRUSKAL–WALLIS TABLE

results = [];

for i = 1:numel(groupNames)
    for j = i+1:numel(groupNames)

        dataA = groupData{i};
        dataB = groupData{j};

        % Skip empty groups (should not happen now)
        if isempty(dataA) || isempty(dataB)
            continue
        end

        combined = [dataA; dataB];
        labels = [repmat("A", numel(dataA), 1); repmat("B", numel(dataB), 1)];

        pKW = kruskalwallis(combined, labels, 'off');

        results = [results; {
            groupNames(i), groupNames(j), pKW
        }];
    end
end

KW_Pairwise_Stiffness = cell2table(results, ...
    'VariableNames', {'GroupA','GroupB','pValue'});

%% STATISTICAL ANALYSIS — SWELLING (percentM)

fprintf("\n=== STATISTICAL ANALYSIS (Swelling: percentM) ===\n");


% Containers for normality results
normalityMatrix_swelling = false(numTimes, numTemps);

% --- NORMALITY TESTING LOOP ---
for i = 1:numTimes
    for j = 1:numTemps

        % Logical index for this group
        idx = (crosslinkTime == uniqueTimes(i)) & ...
              (crosslinkTemp == uniqueTemps(j));

        dataTemp = percentM(idx);   % <-- swelling data

        if isempty(dataTemp)
            continue
        end

        % Shapiro–Wilk normality test
        [pSW, ~] = swft(dataTemp);
        normalityMatrix_swelling(i,j) = pSW{2,7} > 0.05;

    end
end

% --- GLOBAL NORMALITY DECISION ---
allNormal_swelling = all(normalityMatrix_swelling(:));
if allNormal_swelling
    fprintf("\nNormal groups detected (Swelling).\n");
else
    fprintf("\nNon‑normal groups detected (Swelling) > Running KRUSKAL–WALLIS.\n");
end


% ============================================================
%   PAIRWISE KRUSKAL–WALLIS FOR SWELLING (percentM)
%   Using the same 7 groups as stiffness analysis
% ============================================================

fprintf("\n=== Pairwise Kruskal–Wallis Tests (Swelling) ===\n");

% 0 hr (uncrosslinked) — no temperature distinction
S_0h = percentM(inventory.CrosslinkTime == 0);

% 2 hr groups
S_2h_27 = percentM(inventory.CrosslinkTime == 2 & inventory.CrosslinkTemp == 27);
S_2h_37 = percentM(inventory.CrosslinkTime == 2 & inventory.CrosslinkTemp == 37);

% 4 hr groups
S_4h_27 = percentM(inventory.CrosslinkTime == 4 & inventory.CrosslinkTemp == 27);
S_4h_37 = percentM(inventory.CrosslinkTime == 4 & inventory.CrosslinkTemp == 37);

% 6 hr groups
S_6h_27 = percentM(inventory.CrosslinkTime == 6 & inventory.CrosslinkTemp == 27);
S_6h_37 = percentM(inventory.CrosslinkTime == 6 & inventory.CrosslinkTemp == 37);

% REGISTER GROUPS
groupNames_swelling = [
    "0h"
    "2h_27"
    "2h_37"
    "4h_27"
    "4h_37"
    "6h_27"
    "6h_37"
];

groupData_swelling = {
    S_0h
    S_2h_27
    S_2h_37
    S_4h_27
    S_4h_37
    S_6h_27
    S_6h_37
};

% PAIRWISE KRUSKAL–WALLIS TABLE (SWELLING)

results_swelling = [];

for i = 1:numel(groupNames_swelling)
    for j = i+1:numel(groupNames_swelling)

        dataA = groupData_swelling{i};
        dataB = groupData_swelling{j};

        if isempty(dataA) || isempty(dataB)
            continue
        end

        combined = [dataA; dataB];
        labels   = [repmat("A", numel(dataA), 1); repmat("B", numel(dataB), 1)];

        pKW = kruskalwallis(combined, labels, 'off');

        results_swelling = [results_swelling; {
            groupNames_swelling(i), groupNames_swelling(j), pKW
        }];
    end
end

KW_Pairwise_Swelling = cell2table(results_swelling, ...
    'VariableNames', {'GroupA','GroupB','pValue'});


%% PLOT STIFFNESS GROUPS

% Create the boxchart
stiffness = figure;

% Plot the data using boxchart, with crosslink time on the x-axis and effective modulus as the y-axis
% Group by color based on temperature
boxchart(xShiftedBox, modulusData, 'GroupByColor', crosslinkTemp, 'BoxWidth', 1);
colororder(colors);
hold on; 

scatter(1, meanE0, "Marker", "o", "MarkerFaceColor","#000000");
scatter([1.85, 2.85, 3.85], meanE27, "Marker", "o", "MarkerFaceColor","#4DBEEE");
scatter([2.15, 3.15, 4.15], meanE37, "Marker", "o", "MarkerFaceColor","#D95319");

gscatter(xShifted, modulusData, crosslinkTemp, colors, 'o'); 

% Adjust the labels and title
xlabel('Crosslink Time (hr)');
ylabel('Effective Elastic Modulus (kPa)');
title('Effective Elastic Modulus vs Crosslink Time');


% Adjust X-tick labels if needed (e.g., if you want to display them as '2hr', '4hr', '6hr')
xticks(positions);
xticklabels(strcat(string(unique(crosslinkTime)), 'hr'));

H27 = sigstar({[1, 1.85]}, 0.0260, 1);
% sigstar({[1.875, 2.875]}, 0.8182);
% sigstar({[2.875, 3.875]}, 0.3939);
set(H27,'color','#4DBEEE');

H37 = sigstar({[1, 2.15]}, 0.0260, 1);
% sigstar({[2.125, 3.125]}, 0.3939);
% sigstar({[3.125, 4.125]}, 0.3939);
set(H37,'color','#D95319');

%H2737 = sigstar({[1.85, 2.15], [2.875, 3.125], [3.875, 4.125]}, [0.5887, 0.8182, 0.9372], 1);
%set(H2737,'color','k');

% Add a legend for the temperature groups
legend('baseline', '27 C', '37 C', 'Location', 'southeast');

% Customize the appearance if desired (optional)
set(gca, 'FontSize', 12); % Adjust font size for better readability
box on; % Draw box around the plot

pbaspect([3, 3, 1]);

hold off;
%saveas(stiffness, root + "/FIGURES/Effective Elastic Modulus vs Crosslink Time.svg", 'svg');


%% PLOT NET MASS CONTRACTION GROUPS

% Create the boxchart
mass = figure;
hold on;

% Plot the data using boxchart, with crosslink time on the x-axis and effective modulus as the y-axis
% Group by color based on temperature
boxchart(xShiftedBox, percentM, 'GroupByColor', crosslinkTemp);
colororder(colors);

scatter([1.85, 2.85, 3.85], meanPM27, "Marker", "o", "MarkerFaceColor","#4DBEEE");
scatter([2.15, 3.15, 4.15], meanPM37, "Marker", "o", "MarkerFaceColor","#D95319");

gscatter(xShifted, percentM, crosslinkTemp, colors, 'o'); 

% Adjust the labels and title
xlabel('Crosslink Time (hr)');
ylabel('Negative Swelling (% Change in Mass)');
title('Mass Swelling vs Crosslink Time');

% Adjust X-tick labels if needed (e.g., if you want to display them as '2hr', '4hr', '6hr')
xticks(positions);
xticklabels(strcat(string(unique(crosslinkTime)), 'hr'));

HM27 = sigstar({[1, 1.85]}, 0.0022, 1);
% sigstar({[1.875, 2.875]}, 0.2403);
% sigstar({[2.875, 3.875]}, 0.3095);
set(HM27,'color','#4DBEEE');

HM37 = sigstar({[1, 2.15]}, 0.0022, 1);
% sigstar({[2.125, 3.125]}, 0.3095);
% sigstar({[3.125, 4.125]}, 0.3095);
set(HM37,'color','#D95319');

H2737 = sigstar({[2.85, 3.15], [3.85, 4.15]}, [0.0022, 0.0152], 1);
% sigstar({[1.875, 2.125]}, 0.3939);
set(H2737,'color','k');

% Add a legend for the temperature groups
legend('baseline', '27 C', '37 C', 'Location', 'southeast');

% Customize the appearance if desired (optional)
set(gca, 'FontSize', 12); % Adjust font size for better readability
box on; % Draw box around the plot

pbaspect([3, 3, 1]);
currentFolder = pwd;

%saveas(mass, root + "/FIGURES/Mass Swelling vs Crosslink Time Group.svg", 'svg');


%% FITTED MODELS DATA
% 
% % select columns to be plotted
% results = table(inventory.SampleID, inventory.CrosslinkTime, inventory.CrosslinkTemp, inventory.EffectiveElasticModulus, percentM, 'VariableNames', {'SampleID', 'CrosslinkTime', 'CrosslinkTemp', 'ModulusData', 'PercentSM'});
% 
% % Sort data by crosslinkTemp
% sortedResults = sortrows(results, 'CrosslinkTemp');
% 
% % Separate data by temperature
% results27 = sortedResults(sortedResults.CrosslinkTemp == 27, :);
% results37  = sortedResults(sortedResults.CrosslinkTemp == 37, :);
% 
% % Include entries where CrosslinkTime is 0 in both datasets
% results27 = [results27; results(results.CrosslinkTime == 0, :)];
% results37 = [results37; results(results.CrosslinkTime == 0, :)];
% 
% % Remove duplicate rows
% results27 = unique(results27, 'rows');
% results37 = unique(results37, 'rows');
% 
% %% PLOT FITTED MODELS
% 
% ft = fittype('c - a*exp(-b*x)');
% opts = fitoptions(ft);
% opts.Lower = [-Inf -Inf -Inf];
% opts.Upper = [Inf Inf Inf];
% opts.StartPoint = [0.01 0.5 0.2];
% 
% [fitresult27E, gof27E] = fit(results27.CrosslinkTime, results27.ModulusData, ft, opts);
% 
% % [fitresult37E, gof37E] = fit(results37.CrosslinkTime, results37.ModulusData, 'exp2', 'StartPoint', [-0.06 1 0.1], [-Inf, -Inf, -Inf], [Inf, fitresult27E.b, Inf]);
% % 
% % [fitresult37SM, gof37SM] = createFit(results37.CrosslinkTime, results37.PercentSM, [-10 1 10]);
% opts.StartPoint = [10 1 10];
% [fitresult27SM, gof27SM] = fit(results27.CrosslinkTime, results27.PercentSM, ft, opts);
% % 
% time = linspace(0, 6, 100);
% fitted27E = fitresult27E(time);
% % fitted37E = fitresult37E(time);
% fitted27SM = fitresult27SM(time);
% % fitted37SM = fitresult37SM(time);
% 
% 
% 
% modulusFits = figure;
% hold on;
% scatter(results27.CrosslinkTime, results27.ModulusData, 'filled');
% plot(time, fitted27E, 'Color', 'k');
% xlabel('Crosslink Time (hr)');
% ylabel('Effective Elastic Modulus (Pa)');
% title('Fitted Effective Elastic Modulus vs Crosslink Time');
% txt = "y(x) = 0.179 - 0.062e^{-0.648x}" + newline + "R^2 = " + gof27E.rsquare;
% legend('27C', txt, 'Location', 'southeast');
% 
% saveas(modulusFits, root + "/FIGURES/Fitted Effective Elastic Modulus vs Crosslink Time.svg", 'svg');
% 
% swellingFits = figure;
% hold on;
% scatter(results27.CrosslinkTime, results27.PercentSM, 'filled');
% plot(time, fitted27SM, 'Color', 'k');
% xlabel('Crosslink Time (hr)');
% ylabel('Negative Swelling (%)');
% title('Negative Swelling vs Crosslink Time');
% txt = "y(x) = 13.61 - 12.86e^{-13.61x}" + newline + "R^2 = " + gof27SM.rsquare;
% legend('27C', txt, 'Location', 'southeast');
% 
% saveas(swellingFits, root + "/FIGURES/Fitted Negative Swelling vs Crosslink Time.svg", 'svg');
% 
% 

%% PLOT MONOMOER RATIOS (literature) 
% 
% % Data for the bar chart
% species = {'SAMPLE 1', 'SAMPLE 2', 'SAMPLE 3', 'Sargassum asperifolium', 'Sargassum duplicatum', 'Sargassum ilicifolium',  'Sargassum patens', ...
%            'Sargassum olgiocystum', 'Sargassum vulgare', 'Turbinaria ornata', 'Padina gymnospora', 'Cystoseira compressa', 'Laminaria digitata'};
% mg_ratio = [0.55, 0.60, 0.91, 0.69, 0.86, 0.25, 1.59, 0.49, 0.71, 1.25, 1.5, 0.77, 1.12];
% 
% % Set species names with some in italics and some in bold
% italic_species = strcat('\it', species);
% bold_species = strcat('\bf', species);
% formatted_species = italic_species; % Default to italics
% 
% % Normalize mg_ratio for color mapping
% norm_mg_ratio = (mg_ratio - min(mg_ratio)) / (max(mg_ratio) - min(mg_ratio));
% 
% % Create the figure
% figure('Position', [100, 100, 400, 500]); % Narrow and tall figure
% 
% % Create the scatter plot
% monomerComposition = scatter(mg_ratio, 1:numel(species), 100, norm_mg_ratio, 'filled');
% colormap([linspace(0, 1, 256)', linspace(0, 0, 256)', linspace(1, 0, 256)']); % blue to red colormap
% colorbar('Ticks', [], 'TickLabels', [], 'Location', 'southoutside', 'Orientation', 'horizontal'); % Remove ticks and labels from colorbar
% caxis([0 1]);
% 
% % Set species names with some in italics and some in bold
% italic_species = strcat('\it', species);
% bold_species = strcat('\bf', species);
% formatted_species = italic_species; % Default to italics
% 
% % Specify which species to bold
% bold_indices = [1, 2, 3]; % Example indices to be bolded
% formatted_species(bold_indices) = bold_species(bold_indices);
% 
% % Set the y-axis labels
% set(gca, 'YTick', 1:numel(species), 'YTickLabel', formatted_species);
% set(gca, 'YDir', 'reverse'); % Invert y-axis to match the original orientation
% 
% % Add title and axis labels
% title('Monomer Composition of Known Alginates');
% xlabel('M/G Ratio');
% ylabel('Species');
% 
% xlim([0 2]);
% ytickangle(45);
% 
% saveas(monomerComposition, root + "/FIGURES/Monomer Composition of Alginates.svg", 'svg');
% 
%% PLOT MONOMER RATIOS (add samples from HNMR)

% Select data points from indices 4 to 13 (MATLAB indexing starts at 1)
selected_indices = 4:13;
selected_species = species(selected_indices);
selected_mg_ratio = mg_ratio(selected_indices);

% Normalize mg_ratio for color mapping
norm_mg_ratio = (selected_mg_ratio - min(selected_mg_ratio)) / (max(selected_mg_ratio) - min(selected_mg_ratio));

% Create the figure
figure('Position', [100, 100, 400, 500]); % Narrow and tall figure

% Create the scatter plot
monomerComposition2 = scatter(selected_mg_ratio, 1:numel(selected_species), 100, norm_mg_ratio, 'filled');
colormap([linspace(0, 1, 256)', linspace(0, 0, 256)', linspace(1, 0, 256)']); % blue to red colormap
colorbar('Ticks', [], 'TickLabels', [], 'Location', 'southoutside', 'Orientation', 'horizontal'); % Remove ticks and labels from colorbar
caxis([0 1]);

% Set species names with some in italics and some in bold
italic_species = strcat('\it', selected_species);
bold_species = strcat('\bf', selected_species);
formatted_species = italic_species; % Default to italics

% Set the y-axis labels
set(gca, 'YTick', 1:numel(selected_species), 'YTickLabel', formatted_species);
set(gca, 'YDir', 'reverse'); % Invert y-axis to match the original orientation

% Add title and axis labels
title('Monomer Composition of Known Alginates');
xlabel('M/G Ratio');
ylabel('Species');

% Set x-axis limits
xlim([0 2]);
ytickangle(45);
ylim([0 numel(selected_species) + 1]);

% saveas(monomerComposition2, root + "/FIGURES/Monomer Composition of Alginates (no samples).svg", 'svg');




%% PLOT INITIAL GEOMETRY 


% Create the boxchart
initialGeom = figure;
hold on;

% Plot the data using boxchart, with crosslink time on the x-axis and effective modulus as the y-axis
% Group by color based on temperature

boxchart(xShiftedBox, initialA, 'GroupByColor', crosslinkTemp, 'BoxWidth', 1);
colororder(colors);
hold on; 


scatter(1, meanA0, "Marker", "o", "MarkerFaceColor","#000000");
scatter([1.85, 2.85, 3.85], meanA27, "Marker", "o", "MarkerFaceColor","#4DBEEE");
scatter([2.15, 3.15, 4.15], meanA37, "Marker", "o", "MarkerFaceColor","#D95319");

gscatter(xShifted, initialA, crosslinkTemp, colors, 'o'); 


% Adjust the labels and title
xlabel('Crosslink Time Group');
ylabel('Cross Sectional Area (mm^2)');
title('Initial Cross Sectional Area');

% Adjust X-tick labels if needed (e.g., if you want to display them as '2hr', '4hr', '6hr')
xticks(positions);
xticklabels(strcat(string(unique(crosslinkTime)), 'hr'));

% Customize the appearance if desired (optional)
set(gca, 'FontSize', 12); % Adjust font size for better readability
box on; % Draw box around the plot
pbaspect([3, 3, 1]);

% Add a legend for the temperature groups
legend('baseline', '27 C', '37 C', 'mean', 'Location', 'northwest');
hold off;
saveas(initialGeom, root + "/FIGURES/Initial Geometry.svg", 'svg');

%% PLOT INITIAL MASS 


% Create the boxchart
initialMass = figure;
hold on;

% Plot the data using boxchart, with crosslink time on the x-axis and effective modulus as the y-axis
% Group by color based on temperature

boxchart(xShiftedBox, initialM, 'GroupByColor', crosslinkTemp, 'BoxWidth', 1);
colororder(colors);
hold on; 


scatter(1, meanM0, "Marker", "o", "MarkerFaceColor","#000000");
scatter([1.85, 2.85, 3.85], meanM27, "Marker", "o", "MarkerFaceColor","#4DBEEE");
scatter([2.15, 3.15, 4.15], meanM37, "Marker", "o", "MarkerFaceColor","#D95319");

gscatter(xShifted, initialM, crosslinkTemp, colors, 'o'); 


% Adjust the labels and title
xlabel('Crosslink Time Group');
ylabel('Mass (g)');
title('Initial Mass');

% Adjust X-tick labels if needed (e.g., if you want to display them as '2hr', '4hr', '6hr')
xticks(positions);
xticklabels(strcat(string(unique(crosslinkTime)), 'hr'));

% Customize the appearance if desired (optional)
set(gca, 'FontSize', 12); % Adjust font size for better readability
box on; % Draw box around the plot
pbaspect([3, 3, 1]);

% Add a legend for the temperature groups
legend('baseline', '27 C', '37 C', 'mean', 'Location', 'northwest');
hold off;
saveas(initialMass, root + "/FIGURES/Initial Mass.svg", 'svg');

