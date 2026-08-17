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

initialM = inventory.initialMass_g_;
finalM = inventory.finalMass_g_;

meanModulus = groupsummary(modulusData, [crosslinkTime, crosslinkTemp],'mean');
mean0 = [meanModulus(1)];
mean27 = [meanModulus(2); meanModulus(4); meanModulus(6)];
mean37 = [meanModulus(3); meanModulus(5); meanModulus(7)];

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
temperatureOffset = [0, -0.125, 0.125]; % Offset for '27°C' and '37°C' to avoid overlap
temperatureOffsetBox = [0.333, -0.125, -0.208];
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
boxchart(xShiftedBox, modulusData, 'GroupByColor', crosslinkTemp);
colororder(colors);
hold on; 

scatter(1, mean0, "Marker", "o", "MarkerFaceColor","#000000");
scatter([1.875, 2.875, 3.875], mean27, "Marker", "o", "MarkerFaceColor","#4DBEEE");
scatter([2.125, 3.125, 4.125], mean37, "Marker", "o", "MarkerFaceColor","#D95319");

gscatter(xShifted, modulusData, crosslinkTemp, [], 'o'); 

% Adjust the labels and title
xlabel('Crosslink Time (hr)');
ylabel('Effective Elastic Modulus (kPa)');
title('Effective Elastic Modulus vs Crosslink Time');


% Adjust X-tick labels if needed (e.g., if you want to display them as '2hr', '4hr', '6hr')
xticks(positions);
xticklabels(strcat(string(unique(crosslinkTime)), 'hr'));

H27 = sigstar({[1, 1.875]}, 0.0260, 1);
% sigstar({[1.875, 2.875]}, 0.8182);
% sigstar({[2.875, 3.875]}, 0.3939);
set(H27,'color','#4DBEEE');

H37 = sigstar({[1, 2.125]}, 0.0260, 1);
% sigstar({[2.125, 3.125]}, 0.3939);
% sigstar({[3.125, 4.125]}, 0.3939);
set(H37,'color','#D95319');

%H2737 = sigstar({[1.875, 2.125], [2.875, 3.125], [3.875, 4.125]}, [0.5887, 0.8182, 0.9372], 1);
%set(H2737,'color','k');

% Add a legend for the temperature groups
legend('uncrosslinked', '27 C', '37 C', 'Location', 'southeast');

% Customize the appearance if desired (optional)
set(gca, 'FontSize', 12); % Adjust font size for better readability
box on; % Draw box around the plot

%pbaspect([3, 3, 1]);

hold off;
%saveas(stiffness, root + "MTS Experiments/Ring/FIGURES/Effective Elastic Modulus vs Crosslink Time.svg", 'svg');


%% PLOT NET SWELLING GROUPS

% Create the boxchart
mass = figure;
hold on;

% Plot the data using boxchart, with crosslink time on the x-axis and effective modulus as the y-axis
% Group by color based on temperature
boxchart(xShiftedBox, percentM, 'GroupByColor', crosslinkTemp);
colororder(colors);

scatter([1.875, 2.875, 3.875], meanPM27, "Marker", "o", "MarkerFaceColor","#4DBEEE");
scatter([2.125, 3.125, 4.125], meanPM37, "Marker", "o", "MarkerFaceColor","#D95319");

gscatter(xShifted, percentM, crosslinkTemp, [], 'o'); 

% Adjust the labels and title
xlabel('Crosslink Time (hr)');
ylabel('Negative Swelling (% Change in Mass)');
title('Mass Swelling vs Crosslink Time');

% Adjust X-tick labels if needed (e.g., if you want to display them as '2hr', '4hr', '6hr')
xticks(positions);
xticklabels(strcat(string(unique(crosslinkTime)), 'hr'));

HM27 = sigstar({[1, 1.875]}, 0.0022, 1);
% sigstar({[1.875, 2.875]}, 0.2403);
% sigstar({[2.875, 3.875]}, 0.3095);
set(HM27,'color','#4DBEEE');

HM37 = sigstar({[1, 2.125]}, 0.0022, 1);
% sigstar({[2.125, 3.125]}, 0.3095);
% sigstar({[3.125, 4.125]}, 0.3095);
set(HM37,'color','#D95319');

H2737 = sigstar({[2.875, 3.125], [3.875, 4.125]}, [0.0022, 0.0152], 1);
% sigstar({[1.875, 2.125]}, 0.3939);
set(H2737,'color','k');

% Add a legend for the temperature groups
legend('uncrosslinked', '27 C', '37 C', 'Location', 'southeast');

% Customize the appearance if desired (optional)
set(gca, 'FontSize', 12); % Adjust font size for better readability
box on; % Draw box around the plot

%pbaspect([3, 3, 1]);
currentFolder = pwd;

saveas(mass, root + "/Experimental Data/FIGURES/Mass Swelling vs Crosslink Time Group.svg", 'svg');

%% FITTED MODELS DATA

% select columns to be plotted
results = table(inventory.SampleID, inventory.CrosslinkTime, inventory.CrosslinkTemp, inventory.EffectiveElasticModulus, percentM, 'VariableNames', {'SampleID', 'CrosslinkTime', 'CrosslinkTemp', 'ModulusData', 'PercentSM'});

% Sort data by crosslinkTemp
sortedResults = sortrows(results, 'CrosslinkTemp');

% Separate data by temperature
results27 = sortedResults(sortedResults.CrosslinkTemp == 27, :);
results37  = sortedResults(sortedResults.CrosslinkTemp == 37, :);

% Include entries where CrosslinkTime is 0 in both datasets
results27 = [results27; results(results.CrosslinkTime == 0, :)];
results37 = [results37; results(results.CrosslinkTime == 0, :)];

% Remove duplicate rows
results27 = unique(results27, 'rows');
results37 = unique(results37, 'rows');

%% PLOT FITTED MODELS

ft = fittype('c - a*exp(-b*x)');
opts = fitoptions(ft);
opts.Lower = [-Inf -Inf -Inf];
opts.Upper = [Inf Inf Inf];
opts.StartPoint = [0.01 0.5 0.2];

[fitresult27E, gof27E] = fit(results27.CrosslinkTime, results27.ModulusData, ft, opts);

% [fitresult37E, gof37E] = fit(results37.CrosslinkTime, results37.ModulusData, 'exp2', 'StartPoint', [-0.06 1 0.1], [-Inf, -Inf, -Inf], [Inf, fitresult27E.b, Inf]);
% 
% [fitresult37SM, gof37SM] = createFit(results37.CrosslinkTime, results37.PercentSM, [-10 1 10]);
opts.StartPoint = [10 1 10];
[fitresult27SM, gof27SM] = fit(results27.CrosslinkTime, results27.PercentSM, ft, opts);
% 
time = linspace(0, 6, 100);
fitted27E = fitresult27E(time);
% fitted37E = fitresult37E(time);
fitted27SM = fitresult27SM(time);
% fitted37SM = fitresult37SM(time);



modulusFits = figure;
hold on;
scatter(results27.CrosslinkTime, results27.ModulusData, 'filled');
plot(time, fitted27E, 'Color', 'k');
xlabel('Crosslink Time (hr)');
ylabel('Effective Elastic Modulus (Pa)');
title('Fitted Effective Elastic Modulus vs Crosslink Time');
txt = "y(x) = 0.179 - 0.062e^{-0.648x}" + newline + "R^2 = " + gof27E.rsquare;
legend('27C', txt, 'Location', 'southeast');

saveas(modulusFits, root + "/Experimental Data/FIGURES/Fitted Effective Elastic Modulus vs Crosslink Time.svg", 'svg');

swellingFits = figure;
hold on;
scatter(results27.CrosslinkTime, results27.PercentSM, 'filled');
plot(time, fitted27SM, 'Color', 'k');
xlabel('Crosslink Time (hr)');
ylabel('Negative Swelling (%)');
title('Negative Swelling vs Crosslink Time');
txt = "y(x) = 13.61 - 12.86e^{-13.61x}" + newline + "R^2 = " + gof27SM.rsquare;
legend('27C', txt, 'Location', 'southeast');

saveas(swellingFits, root + "/Experimental Data/FIGURES/Fitted Negative Swelling vs Crosslink Time.svg", 'svg');


%% PROCESS DYNAMIC SWELLING DATA 
A_file = "Areas Final new.csv";
A = readmatrix(root+A_file);

skipped = A(:, 1) == 0;
A(skipped, :) = [];


% for i = 1:43
%     A(i, :) = sort(A(i, :), 'descend');
% end

Anorm = A(:, :)./A(:, 1) * 100;

Anorm(Anorm==0) = nan;
writematrix(Anorm, 'Anorm.csv');

sortedResults = sortrows(inventory, 'CrosslinkTime');

S2_indx = (inventory.CrosslinkTime == 2 & inventory.CrosslinkTemp == 27);
S4_indx = (inventory.CrosslinkTime == 4 & inventory.CrosslinkTemp == 27);
S6_indx = (inventory.CrosslinkTime == 6 & inventory.CrosslinkTemp == 27);
T2_indx = (inventory.CrosslinkTime == 2 & inventory.CrosslinkTemp == 37);
T4_indx = (inventory.CrosslinkTime == 4 & inventory.CrosslinkTemp == 37);
T6_indx = (inventory.CrosslinkTime == 6 & inventory.CrosslinkTemp == 37);


% Define time variable
time2 = [0, 15, 30, 45, 60, 75, 90, 105, 120];
time4 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240];
time6 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 360];

rep = 6; % replicates

S2 = [Anorm(S2_indx, 1:9)];
fitS2 = zeros(length(S2(:, 1)), 9);
S4 = [Anorm(S4_indx, 1:17)];
fitS4 = zeros(length(S4(:, 1)), 17);
S6 = [Anorm(S6_indx, 1:25)];
fitS6 = zeros(length(S6(:, 1)), 25);
T2 = [Anorm(T2_indx, 1:9)];
fitT2 = zeros(length(T2(:, 1)), 9);
T4 = [Anorm(T4_indx, 1:17)];
fitT4 = zeros(length(T4(:, 1)), 17);
T6 = [Anorm(T6_indx, 1:25)];
fitT6 = zeros(length(T6(:, 1)), 25);

y0 = 100; 
% Set up fittype and options.
ft = fittype( 'a1*exp(-b1*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
opts.StartPoint = [0.349983765984809 0.196595250431208 90];

for i = 1:rep

    fS2 = fit(time2', S2(i, :)', ft, opts);
    fS4 = fit(time4', S4(i, :)', ft, opts);
    fS6 = fit(time6', S6(i, :)', ft, opts);
    fT2 = fit(time2', T2(i, :)', ft, opts);
    fT4 = fit(time4', T4(i, :)', ft, opts);
    fT6 = fit(time6', T6(i, :)', ft, opts);


    % Adjust each fitted model to start at y0 by adding a constant shift
    % Calculate the adjustment at x = 0 (initial values)
    adjustmentS2 = y0 - fS2(0);
    adjustmentS4 = y0 - fS4(0);
    adjustmentS6 = y0 - fS6(0);
    adjustmentT2 = y0 - fT2(0);
    adjustmentT4 = y0 - fT4(0);
    adjustmentT6 = y0 - fT6(0);

    % Store the adjusted fitted values
    fitS2(i, :) = fS2(time2')' + adjustmentS2;
    fitS4(i, :) = fS4(time4')' + adjustmentS4;
    fitS6(i, :) = fS6(time6')' + adjustmentS6;
    fitT2(i, :) = fT2(time2')' + adjustmentT2;
    fitT4(i, :) = fT4(time4')' + adjustmentT4;
    fitT6(i, :) = fT6(time6')' + adjustmentT6;


end

mS2 = mean(S2);
mfitS2 = zeros(1, 9);
mS4 = mean(S4);
mfitS4 = zeros(1, 17);
mS6 = mean(S6);
mfitS6 = zeros(1, 25);
mT2 = mean(T2);
mfitT2 = zeros(1, 9);
mT4 = mean(T4);
mfitT4 = zeros(1, 17);
% for T6, ignore outlier
%tempT6 = [T6(1, :); T6(3, :)];
mT6 = mean(T6);
mfitT6 = zeros(1, 25);


    mfS2 = fit(time2', mS2', ft, opts);
    mfS4 = fit(time4', mS4', ft, opts);
    mfS6 = fit(time6', mS6', ft, opts);
    mfT2 = fit(time2', mT2', ft, opts);
    mfT4 = fit(time4', mT4', ft, opts);
    mfT6 = fit(time6', mT6', ft, opts);

    % Adjust each fitted model to start at y0 by adding a constant shift
    % Calculate the adjustment at x = 0 (initial values)
    madjustmentS2 = y0 - mfS2(0);
    madjustmentS4 = y0 - mfS4(0);
    madjustmentS6 = y0 - mfS6(0);
    madjustmentT2 = y0 - mfT2(0);
    madjustmentT4 = y0 - mfT4(0);
    madjustmentT6 = y0 - mfT6(0);

    % Store the adjusted fitted values
    mfitS2 = mfS2(time2')' + madjustmentS2;
    mfitS4 = mfS4(time4')' + madjustmentS4;
    mfitS6 = mfS6(time6')' + madjustmentS6;
    mfitT2 = mfT2(time2')' + madjustmentT2;
    mfitT4 = mfT4(time4')' + madjustmentT4;
    mfitT6 = mfT6(time6')' + madjustmentT6;

stdS2 = std(S2);
stdS4 = std(S4);
stdS6 = std(S6);
stdT2 = std(T2);
stdT4 = std(T4);
stdT6 = std(T6);

%% PLOT DYNAMIC SWELLING (all samples, one plot) 

areaOverTime = figure; 
hold on;

% Define colors and line styles
colors = {'r', 'g', 'b'}; % Colors for 2, 4, and 6 time groups
lineStyles = {'-', '--'}; % Line styles for S and T groups

% Plot each group with different colors and styles, hiding extra entries in the legend
for i = 1:rep
    % S groups
    plot(time2, fitS2(i, :), 'Color', colors{1}, 'LineStyle', lineStyles{1}, 'HandleVisibility', 'off');
    plot(time4, fitS4(i, :), 'Color', colors{2}, 'LineStyle', lineStyles{1}, 'HandleVisibility', 'off');
    plot(time6, fitS6(i, :), 'Color', colors{3}, 'LineStyle', lineStyles{1}, 'HandleVisibility', 'off');

    % T groups
    plot(time2, fitT2(i, :), 'Color', colors{1}, 'LineStyle', lineStyles{2}, 'HandleVisibility', 'off');
    plot(time4, fitT4(i, :), 'Color', colors{2}, 'LineStyle', lineStyles{2}, 'HandleVisibility', 'off');
    plot(time6, fitT6(i, :), 'Color', colors{3}, 'LineStyle', lineStyles{2}, 'HandleVisibility', 'off');
end

% Add one line per legend entry manually
plot(NaN, NaN, 'Color', 'r', 'LineStyle', '-', 'DisplayName', '2hr @ 27C');
plot(NaN, NaN, 'Color', 'g', 'LineStyle', '-', 'DisplayName', '4hr @ 27C');
plot(NaN, NaN, 'Color', 'b', 'LineStyle', '-', 'DisplayName', '6hr @ 27C');
plot(NaN, NaN, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '2hr @ 37C');
plot(NaN, NaN, 'Color', 'g', 'LineStyle', '--', 'DisplayName', '4hr @ 37C');
plot(NaN, NaN, 'Color', 'b', 'LineStyle', '--', 'DisplayName', '6hr @ 37C');

% Customize the plot
xlabel('Time (min)');
ylabel('Normalized Area (mm^2)');
title('Exponential Fits for Normalized Area over Time');
%legend('show', 'Location', 'northeast'); % Display the legend with group labels
grid off;
hold off;

% saveas(areaOverTime, root + "/Experimental Data/FIGURES/Normalized Area Over Time (all samples combined).svg", 'svg');

%% PLOT DYNAMIC SWELLING (sample means, one plot) 

areaOverTime = figure; 
hold on;

% Define colors and line styles
colors = {'r', 'g', 'b'}; % Colors for 2, 4, and 6 time groups
lineStyles = {'-', '--'}; % Line styles for S and T groups

% Plot each group with different colors and styles

fill([time2, fliplr(time2)], [mfitS2 + stdS2, fliplr(mfitS2 - stdS2)], colors{1}, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time2, mfitS2, 'Color', colors{1}, 'LineStyle', lineStyles{1}, 'DisplayName', '2hr @ 27C');
fill([time2, fliplr(time2)], [mfitT2 + stdT2, fliplr(mfitT2 - stdT2)], colors{1}, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time2, mfitT2, 'Color', colors{1}, 'LineStyle', lineStyles{2}, 'DisplayName', '2hr @ 37C');

fill([time4, fliplr(time4)], [mfitS4 + stdS4, fliplr(mfitS4 - stdS4)], colors{2}, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time4, mfitS4, 'Color', colors{2}, 'LineStyle', lineStyles{1}, 'DisplayName', '4hr @ 27C');
fill([time4, fliplr(time4)], [mfitT4 + stdT4, fliplr(mfitT4 - stdT4)], colors{2}, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time4, mfitT4, 'Color', colors{2}, 'LineStyle', lineStyles{2}, 'DisplayName', '4hr @ 37C');

fill([time6, fliplr(time6)], [mfitS6 + stdS6, fliplr(mfitS6 - stdS6)], colors{3}, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time6, mfitS6, 'Color', colors{3}, 'LineStyle', lineStyles{1}, 'DisplayName', '6hr @ 27C');
fill([time6, fliplr(time6)], [mfitT6 + stdT6, fliplr(mfitT6 - stdT6)], colors{3}, 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time6, mfitT6, 'Color', colors{3}, 'LineStyle', lineStyles{2}, 'DisplayName', '6hr @ 37C');

% Customize the plot
xlabel('Time (min)');
ylabel('Normalized Area (mm^2)');
title('Exponential Fits for Normalized Area over Time');
legend('show', 'Location', 'northeast'); % Display the legend with group labels
grid off;
hold off;
xlim([0, 365]);
xticks(linspace(0, 360, 25));
ylim([80, 100]);
yticks(linspace(80, 100, 5));

saveas(areaOverTime, root + "/Experimental Data/FIGURES/Normalized Area Over Time Summary.svg", 'svg');

%% PLOT DYNAMIC SWELLING (sample means, different plots)


% Define colors and line styles
colorshex = ["#cc6677", "#44aa99", "#332288"];
colors = {'r', 'g', 'b'}; % Colors for 2, 4, and 6 time groups
lineStyles = {'-', '--'}; % Line styles for S and T groups

% plot all together

meanAreas = figure("Position", [0, 0, 600, 100]);
all = tiledlayout(meanAreas, 1, 6, 'TileSpacing','loose', 'Padding','loose');
ylabel(all, 'Normalized Area (mm^2)', 'FontSize', 12);

annotation(meanAreas, 'textbox', [0 0.2 1 0.05], ...
    'String', 'Time (min)', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'FontWeight', 'normal', ...
    'EdgeColor', 'none');

annotation(meanAreas, 'textbox', [0 0.8 1 0.05], ...
    'String', 'Exponential Fits for Swelling over Time', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 14, 'FontWeight', 'bold', ...
    'EdgeColor', 'none');

%title(all, 'Exponential Fits for Swelling over Time', 'FontSize', 14, 'FontWeight', 'bold');


% 2 hr groups
hr2 = nexttile([1, 1]);
hold on;
fill([time2, fliplr(time2)], [mfitS2 + stdS2, fliplr(mfitS2 - stdS2)], colors{1}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '2hr +/- std');
plot(time2, mfitS2, 'Color', colors{1}, 'LineStyle', lineStyles{1}, 'DisplayName', '2hr @ 27C');
fill([time2, fliplr(time2)], [mfitT2 + stdT2, fliplr(mfitT2 - stdT2)], colors{1}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time2, mfitT2, 'Color', colors{1}, 'LineStyle', lineStyles{2}, 'DisplayName', '2hr @ 37C');

legend('Location', 'southwest');

title('2 Hr Sample Group Means', 'FontWeight','normal');
xlim([0, 125]);
xticks(linspace(0, 120, 9));
ylim([80, 100]);
yticks(linspace(80, 100, 5));

pbaspect([1 1 1]);
grid off;
hold off;

% 4 hr groups 
hr4 = nexttile([1, 2]);
hold on;
fill([time4, fliplr(time4)], [mfitS4 + stdS4, fliplr(mfitS4 - stdS4)], colors{2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '4hr +/- std');
plot(time4, mfitS4, 'Color', colors{2}, 'LineStyle', lineStyles{1}, 'DisplayName', '4hr @ 27C');
fill([time4, fliplr(time4)], [mfitT4 + stdT4, fliplr(mfitT4 - stdT4)], colors{2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time4, mfitT4, 'Color', colors{2}, 'LineStyle', lineStyles{2}, 'DisplayName', '4hr @ 37C');

legend('Location', 'southwest');

title('4 Hr Sample Group Means', 'FontWeight','normal');
xlim([0, 245]);
xticks(linspace(0, 240, 17));
ylim([80, 100]);
yticks(linspace(80, 100, 5));

pbaspect([2 1 1]);


grid off;
hold off;

% 6 hr groups
hr6 = nexttile([1, 3]);
hold on;
fill([time6, fliplr(time6)], [mfitS6 + stdS6, fliplr(mfitS6 - stdS6)], colors{3}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '6hr +/- std');
plot(time6, mfitS6, 'Color', colors{3}, 'LineStyle', lineStyles{1}, 'DisplayName', '6hr @ 27C');
fill([time6, fliplr(time6)], [mfitT6 + stdT6, fliplr(mfitT6 - stdT6)], colors{3}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time6, mfitT6, 'Color', colors{3}, 'LineStyle', lineStyles{2}, 'DisplayName', '6hr @ 37C');

legend('Location', 'southwest');

title('6 Hr Sample Group Means', 'FontWeight','normal');
xlim([0, 365]);
xticks(linspace(0, 360, 25));
ylim([80, 100]);
yticks(linspace(80, 100, 5));

pbaspect([3 1 1]);

grid off;
hold off;

height = 5; % in inches

set(hr2, 'Units', 'inches', 'Position', [1, 1, height * 1.5, height*1.5]);
set(hr4, 'Units', 'inches', 'Position', [1, 1, height * 2, height]);
set(hr6, 'Units', 'inches', 'Position', [1, 1, height * 3, height]);

% Adjust figure size for readability
set(meanAreas, 'Units', 'normalized', 'Position', [0.1 0.1 0.95 0.5]);

saveas(meanAreas, root + "/Experimental Data/FIGURES/Normalized Area over Time (Sample Means).svg", 'svg');

% plot individually 
% 2hr
hr2 = figure;
hold on;

fill([time2, fliplr(time2)], [mfitS2 + stdS2, fliplr(mfitS2 - stdS2)], colors{1}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '2hr +/- std');
plot(time2, mfitS2, 'Color', colors{1}, 'LineStyle', lineStyles{1}, 'DisplayName', '2hr @ 27C');
fill([time2, fliplr(time2)], [mfitT2 + stdT2, fliplr(mfitT2 - stdT2)], colors{1}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time2, mfitT2, 'Color', colors{1}, 'LineStyle', lineStyles{2}, 'DisplayName', '2hr @ 37C');

legend('Location', 'southwest');

title('2 Hr Sample Group Means', 'FontWeight','normal');
xlabel('Time (min)');
ylabel('Normalized Area (mm^2)');
xlim([0, 125]);
xticks(linspace(0, 120, 9));
ylim([80, 100]);
yticks(linspace(80, 100, 5));

pbaspect([1 1 1]);

grid off;
hold off;


% 4hr
hr4 = figure;
hold on;
fill([time4, fliplr(time4)], [mfitS4 + stdS4, fliplr(mfitS4 - stdS4)], colors{2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '4hr +/- std');
plot(time4, mfitS4, 'Color', colors{2}, 'LineStyle', lineStyles{1}, 'DisplayName', '4hr @ 27C');
fill([time4, fliplr(time4)], [mfitT4 + stdT4, fliplr(mfitT4 - stdT4)], colors{2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time4, mfitT4, 'Color', colors{2}, 'LineStyle', lineStyles{2}, 'DisplayName', '4hr @ 37C');

legend('Location', 'southwest');

title('4 Hr Sample Group Means', 'FontWeight','normal');
xlabel('Time (min)');
ylabel('Normalized Area (mm^2)');
xlim([0, 245]);
xticks(linspace(0, 240, 17));
ylim([80, 100]);
yticks(linspace(80, 100, 5));

pbaspect([2 1 1]);

grid off;
hold off;


% 6hr 
hr6 = figure;
hold on;
fill([time6, fliplr(time6)], [mfitS6 + stdS6, fliplr(mfitS6 - stdS6)], colors{3}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '6hr +/- std');
plot(time6, mfitS6, 'Color', colors{3}, 'LineStyle', lineStyles{1}, 'DisplayName', '6hr @ 27C');
fill([time6, fliplr(time6)], [mfitT6 + stdT6, fliplr(mfitT6 - stdT6)], colors{3}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time6, mfitT6, 'Color', colors{3}, 'LineStyle', lineStyles{2}, 'DisplayName', '6hr @ 37C');

legend('Location', 'southwest');

title('6 Hr Sample Group Means', 'FontWeight','normal');
xlabel('Time (min)');
ylabel('Normalized Area (mm^2)');
xlim([0, 365]);
xticks(linspace(0, 360, 25));
ylim([80, 100]);
yticks(linspace(80, 100, 5));

pbaspect([3 1 1]);

grid off;
hold off;

% save all 

height = 5; % in inches

set(hr2, 'Units', 'inches', 'Position', [1, 1, height * 1, height]);
set(hr4, 'Units', 'inches', 'Position', [1, 1, height * 2, height]);
set(hr6, 'Units', 'inches', 'Position', [1, 1, height * 3, height]);

saveas(hr2, root + "/Experimental Data/FIGURES/Normalized Area over Time (2Hr Sample Means).svg", 'svg');
saveas(hr4, root + "/Experimental Data/FIGURES/Normalized Area over Time (4Hr Sample Means).svg", 'svg');
saveas(hr6, root + "/Experimental Data/FIGURES/Normalized Area over Time (6Hr Sample Means).svg", 'svg')


%% PLOT DYNAMIC SWELLING (specific groups)

time = time2;
G = S2;
label = '2hr @ 27C';


areas_data_per_group = figure;
hold on;

% Plot each group with different colors and styles, hiding extra entries in the legend
fitG = zeros(rep, length(time));
Rs = zeros(rep, 1);
for i = 1:rep

    [fG, gof] = fit(time', G(i, :)', ft, opts);
    Rs(i) = gof.rsquare;

    fitG(i, :) = fG(time')';

    plot(time, fitG(i, :), 'Color', 'r', 'HandleVisibility', 'off');
    plot(time, G(i, :), 'Color', 'k', 'Marker', 'o', 'LineStyle', '--', 'HandleVisibility', 'off');

end

% Add one line per legend entry manually
plot(NaN, NaN, 'Color', 'k', 'Marker', 'o', 'LineStyle', '--', 'DisplayName', label);
label = string({['mean R^2: ' + string(mean(Rs))]});
plot(NaN, NaN, 'Color', 'r', 'DisplayName', label);


% Customize the plot
xlabel('Time (min)');
ylabel('Normalized Area (mm^2)');
title('Exponential Fits for Normalized Area over Time (by group)');
legend('show', 'Location', 'northeast'); % Display the legend with group labels
%text(label, 'northeast')
grid off;
hold off;

saveas(areas_data_per_group, root + "/Experimental Data/FIGURES/Normalized Area over Time (S2).svg", 'svg');


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
% saveas(monomerComposition, root + "/Experimental Data/FIGURES/Monomer Composition of Alginates.svg", 'svg');
% 
%% PLOT MONOMER RATIOS (add samples from HNMR)
% 
% % Select data points from indices 4 to 13 (MATLAB indexing starts at 1)
% selected_indices = 4:13;
% selected_species = species(selected_indices);
% selected_mg_ratio = mg_ratio(selected_indices);
% 
% % Normalize mg_ratio for color mapping
% norm_mg_ratio = (selected_mg_ratio - min(selected_mg_ratio)) / (max(selected_mg_ratio) - min(selected_mg_ratio));
% 
% % Create the figure
% figure('Position', [100, 100, 400, 500]); % Narrow and tall figure
% 
% % Create the scatter plot
% monomerComposition2 = scatter(selected_mg_ratio, 1:numel(selected_species), 100, norm_mg_ratio, 'filled');
% colormap([linspace(0, 1, 256)', linspace(0, 0, 256)', linspace(1, 0, 256)']); % blue to red colormap
% colorbar('Ticks', [], 'TickLabels', [], 'Location', 'southoutside', 'Orientation', 'horizontal'); % Remove ticks and labels from colorbar
% caxis([0 1]);
% 
% % Set species names with some in italics and some in bold
% italic_species = strcat('\it', selected_species);
% bold_species = strcat('\bf', selected_species);
% formatted_species = italic_species; % Default to italics
% 
% % Set the y-axis labels
% set(gca, 'YTick', 1:numel(selected_species), 'YTickLabel', formatted_species);
% set(gca, 'YDir', 'reverse'); % Invert y-axis to match the original orientation
% 
% % Add title and axis labels
% title('Monomer Composition of Known Alginates');
% xlabel('M/G Ratio');
% ylabel('Species');
% 
% % Set x-axis limits
% xlim([0 2]);
% ytickangle(45);
% ylim([0 numel(selected_species) + 1]);
% 
% saveas(monomerComposition2, root + "/Experimental Data/FIGURES/Monomer Composition of Alginates (no samples).svg", 'svg');

%% PLOT NET GEOMETRIC SWELLING 

% percentG = (initialA - finalA) ./ initialA * 100;
% meanPercentG = groupsummary(percentG, [crosslinkTime, crosslinkTemp],'mean');
% meanPG0 = [meanPercentG(1)];
% meanPG27 = [meanPercentG(2); meanPercentG(4); meanPercentG(6)];
% meanPG37 = [meanPercentG(3); meanPercentG(5); meanPercentG(7)];
% 
% % Create the boxchart
% swelling = figure;
% hold on;
% 
% % Plot the data using boxchart, with crosslink time on the x-axis and effective modulus as the y-axis
% % Group by color based on temperature
% boxchart(xShiftedBox, percentG, 'GroupByColor', crosslinkTemp);
% colororder(colors);
% 
% scatter([1.875, 2.875, 3.875], meanPG27, "Marker", "o", "MarkerFaceColor","#4DBEEE");
% scatter([2.125, 3.125, 4.125], meanPG37, "Marker", "o", "MarkerFaceColor","#D95319");
% 
% gscatter(xShifted, percentG, crosslinkTemp, [], 'o'); 
% 
% % Adjust the labels and title
% xlabel('Crosslink Time (hr)');
% ylabel('Negative Swelling (% Change in Area)');
% title('Geometric Swelling vs Crosslink Time');
% 
% % Adjust X-tick labels if needed (e.g., if you want to display them as '2hr', '4hr', '6hr')
% xticks(positions);
% xticklabels(strcat(string(unique(crosslinkTime)), 'hr'));
% 
% % Add a legend for the temperature groups
% legend('uncrosslinked', '27 C', '37 C', 'Location', 'northwest');
% 
% % Customize the appearance if desired (optional)
% set(gca, 'FontSize', 12); % Adjust font size for better readability
% box on; % Draw box around the plot
% 
% 
% saveas(swelling, root + "MTS Experiments/Ring/FIGURES/Geometric Swelling vs Crosslink Time Group.svg", 'svg');



%% PLOT INITIAL GEOMETRY 

% % Create the boxchart
% initialGeom = figure;
% hold on;
% 
% % Plot the data using boxchart, with crosslink time on the x-axis and effective modulus as the y-axis
% % Group by color based on temperature
% b = boxchart(xShiftedBox, initialA, 'GroupByColor', crosslinkTemp);
% colororder(colors);
% 
% 
% % Adjust the labels and title
% xlabel('Crosslink Time (hr)');
% ylabel('Cross Sectional Area (mm^2)');
% title('Initial Geometry vs Sample Group');
% 
% % Adjust X-tick labels if needed (e.g., if you want to display them as '2hr', '4hr', '6hr')
% xticks(positions);
% xticklabels(strcat(string(unique(crosslinkTime)), 'hr'));
% 
% % Customize the appearance if desired (optional)
% set(gca, 'FontSize', 12); % Adjust font size for better readability
% box on; % Draw box around the plot
% 
% gscatter(xShifted, initialA, crosslinkTemp, [], 'o'); 
% 
% % Add a legend for the temperature groups
% legend('uncrosslinked', '27 C', '37 C', 'Location', 'northwest');
% hold off;
% saveas(initialGeom, root + "/Experimental Data/FIGURES/Initial Geometry.svg", 'svg');

