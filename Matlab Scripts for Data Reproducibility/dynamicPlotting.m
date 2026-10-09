
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


initialOD = inventory.initialOD_px_ ./ inventory.Scale_pxpmm_;
initialID = inventory.initialID_px_ ./ inventory.Scale_pxpmm_;

finalOD = inventory.finalOD_px_ ./ inventory.Scale_pxpmm_;
finalID = inventory.finalID_px_ ./ inventory.Scale_pxpmm_;

initialA = (pi()*(initialOD/2).^2) - (pi()*(initialID/2).^2);
finalA = (pi()*(finalOD/2).^2) - (pi()*(finalID/2).^2);

%% PROCESS DYNAMIC SWELLING DATA 
A_file = "Areas Final new.csv";
A = readmatrix(root+A_file);

skipped = A(:, 1) == 0;
A(skipped, :) = [];

Anorm = A(:, :)./A(:, 1);

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


% Set up fittype and options.
ft = fittype( '(1-c)*exp(-b1*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = [0.196595250431208 .90];

for i = 1:rep

    fS2 = fit(time2', S2(i, :)', ft, opts);
    fS4 = fit(time4', S4(i, :)', ft, opts);
    fS6 = fit(time6', S6(i, :)', ft, opts);
    fT2 = fit(time2', T2(i, :)', ft, opts);
    fT4 = fit(time4', T4(i, :)', ft, opts);
    fT6 = fit(time6', T6(i, :)', ft, opts);

    % Store the fitted values
    fitS2(i, :) = fS2(time2')';
    fitS4(i, :) = fS4(time4')';
    fitS6(i, :) = fS6(time6')';
    fitT2(i, :) = fT2(time2')';
    fitT4(i, :) = fT4(time4')';
    fitT6(i, :) = fT6(time6')';

end
% 
% mS2 = mean(S2); mfitS2 = zeros(1, 9);
% mS4 = mean(S4); mfitS4 = zeros(1, 17);
% mS6 = mean(S6); mfitS6 = zeros(1, 25);
% mT2 = mean(T2); mfitT2 = zeros(1, 9);
% mT4 = mean(T4); mfitT4 = zeros(1, 17);
% mT6 = mean(T6); mfitT6 = zeros(1, 25);
% 
% % Fit the mean experimental curves
% mfS2 = fit(time2', mS2', ft, opts);
% mfS4 = fit(time4', mS4', ft, opts);
% mfS6 = fit(time6', mS6', ft, opts);
% mfT2 = fit(time2', mT2', ft, opts);
% mfT4 = fit(time4', mT4', ft, opts);
% mfT6 = fit(time6', mT6', ft, opts);

% % Directly evaluate the mean fit lines
% mfitS2 = mfS2(time2')';
% mfitS4 = mfS4(time4')';
% mfitS6 = mfS6(time6')';
% mfitT2 = mfT2(time2')';
% mfitT4 = mfT4(time4')';
% mfitT6 = mfT6(time6')';

% Calculate experimental means for each group
mS2 = mean(S2, 1, 'omitnan');
mS4 = mean(S4, 1, 'omitnan');
mS6 = mean(S6, 1, 'omitnan');
mT2 = mean(T2, 1, 'omitnan');
mT4 = mean(T4, 1, 'omitnan');
mT6 = mean(T6, 1, 'omitnan');


% Direct assignment: The trendline centerlines are now the empirical means
mfitS2 = mS2;
mfitS4 = mS4;
mfitS6 = mS6;
mfitT2 = mT2;
mfitT4 = mT4;
mfitT6 = mT6;

% Calculate standard deviations for the shaded variance cloud bounds
stdS2 = std(S2, 0, 1, 'omitnan');
stdS4 = std(S4, 0, 1, 'omitnan');
stdS6 = std(S6, 0, 1, 'omitnan');
stdT2 = std(T2, 0, 1, 'omitnan');
stdT4 = std(T4, 0, 1, 'omitnan');
stdT6 = std(T6, 0, 1, 'omitnan');

%% PLOT DYNAMIC SWELLING (all FITTED samples, one plot) 

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

% saveas(areaOverTime, root + "/FIGURES/Normalized Area Over Time (all samples combined).svg", 'svg');

%% PLOT DYNAMIC SWELLING (sample means, NOT FITTED, one plot) 

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

%saveas(areaOverTime, root + "/FIGURES/Normalized Area Over Time Summary.svg", 'svg');

%% PLOT DYNAMIC SWELLING (sample means, different plots)


% Define colors and line styles
colorshex = {'#5e6e3a', '#e99f78', '#f088f7'};
colors = {'r', 'g', 'b'}; % Colors for 2, 4, and 6 time groups
lineStyles = {'-', '--'}; % Line styles for S and T groups

% =========================================================================
% SECTION 1: MASTER TILED LAYOUT PLOT
% =========================================================================

% Initialize figure window dimensionally first
meanAreas = figure("Units", "normalized", "Position", [0.1, 0.1, 0.85, 0.45]);

% Allocate tiles (1 row, 6 columns layout)
allLayout = tiledlayout(meanAreas, 1, 6, 'TileSpacing', 'loose', 'Padding', 'loose');
ylabel(allLayout, 'Normalized Area (mm^2/mm^2)', 'FontSize', 12);

% Add master text annotations over the tiled layout structure
annotation(meanAreas, 'textbox', [0, 0.02, 1, 0.05], ...
    'String', 'Time (min)', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'FontWeight', 'normal', ...
    'EdgeColor', 'none');

annotation(meanAreas, 'textbox', [0, 0.90, 1, 0.05], ...
    'String', 'Exponential Fits for Swelling over Time', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 14, 'FontWeight', 'bold', ...
    'EdgeColor', 'none');


% --- Tile 1: 2 hr groups (#5e6e3a) ---
tile2hr = nexttile(allLayout, [1, 1]);
hold(tile2hr, 'on');

% FIX: Left the 3rd structural position empty ("") so it skips legacy validation
fill(tile2hr, [time2, fliplr(time2)], [mfitS2 + stdS2, fliplr(mfitS2 - stdS2)], "", ...
     'FaceColor', colorshex{1}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '2hr +/- std');

plot(tile2hr, time2, mfitS2, 'Color', colorshex{1}, 'LineStyle', lineStyles{1}, 'DisplayName', '2hr @ 27C');

% FIX: Repeated empty string "" fix here
fill(tile2hr, [time2, fliplr(time2)], [mfitT2 + stdT2, fliplr(mfitT2 - stdT2)], "", ...
     'FaceColor', colorshex{1}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(tile2hr, time2, mfitT2, 'Color', colorshex{1}, 'LineStyle', lineStyles{2}, 'DisplayName', '2hr @ 37C');

legend(tile2hr, 'Location', 'southwest');
title(tile2hr, '2 Hr Sample Group Means', 'FontWeight','normal');
xlim(tile2hr, [0, 125]);
xticks(tile2hr, linspace(0, 120, 9));
ylim(tile2hr, [0.8, 1]); % Adjusted to center around 1.0 normalization
yticks(tile2hr, linspace(0.8, 1, 5));
pbaspect(tile2hr, [1 1 1]);
grid(tile2hr, 'off');


% --- Tile 2: 4 hr groups (#e99f78) ---
tile4hr = nexttile(allLayout, [1, 2]);
hold(tile4hr, 'on');

fill(tile4hr, [time4, fliplr(time4)], [mfitS4 + stdS4, fliplr(mfitS4 - stdS4)], "", ...
     'FaceColor', colorshex{2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '4hr +/- std');

plot(tile4hr, time4, mfitS4, 'Color', colorshex{2}, 'LineStyle', lineStyles{1}, 'DisplayName', '4hr @ 27C');

fill(tile4hr, [time4, fliplr(time4)], [mfitT4 + stdT4, fliplr(mfitT4 - stdT4)], "", ...
     'FaceColor', colorshex{2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(tile4hr, time4, mfitT4, 'Color', colorshex{2}, 'LineStyle', lineStyles{2}, 'DisplayName', '4hr @ 37C');

legend(tile4hr, 'Location', 'southwest');
title(tile4hr, '4 Hr Sample Group Means', 'FontWeight','normal');
xlim(tile4hr, [0, 245]);
xticks(tile4hr, linspace(0, 240, 17));
ylim(tile4hr, [0.8, 1]); 
yticks(tile4hr, linspace(0.8, 1, 5));
pbaspect(tile4hr, [2 1 1]);
grid(tile4hr, 'off');


% --- Tile 3: 6 hr groups (#f088f7) ---
tile6hr = nexttile(allLayout, [1, 3]);
hold(tile6hr, 'on');

fill(tile6hr, [time6, fliplr(time6)], [mfitS6 + stdS6, fliplr(mfitS6 - stdS6)], "", ...
     'FaceColor', colorshex{3}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '6hr +/- std');

plot(tile6hr, time6, mfitS6, 'Color', colorshex{3}, 'LineStyle', lineStyles{1}, 'DisplayName', '6hr @ 27C');

fill(tile6hr, [time6, fliplr(time6)], [mfitT6 + stdT6, fliplr(mfitT6 - stdT6)], "", ...
     'FaceColor', colorshex{3}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');

plot(tile6hr, time6, mfitT6, 'Color', colorshex{3}, 'LineStyle', lineStyles{2}, 'DisplayName', '6hr @ 37C');

legend(tile6hr, 'Location', 'southwest');
title(tile6hr, '6 Hr Sample Group Means', 'FontWeight','normal');
xlim(tile6hr, [0, 365]);
xticks(tile6hr, linspace(0, 360, 25));
ylim(tile6hr, [0.8, 1]);
yticks(tile6hr, linspace(0.8, 1, 5));
pbaspect(tile6hr, [3 1 1]);
grid(tile6hr, 'off');

saveas(meanAreas, root + "/FIGURES/Normalized Area over Time (Sample Means).svg", 'svg');


% =========================================================================
% SECTION 2: STANDALONE INDIVIDUAL FIGURES WITH STRICT PHYSICAL DIMENSIONS
% =========================================================================
height = 5; % baseline height dimension in inches

% --- Standalone 2hr Figure ---
fig2 = figure('Units', 'inches', 'Position', [1, 1, height * 1, height]);
hold on;
fill([time2, fliplr(time2)], [mfitS2 + stdS2, fliplr(mfitS2 - stdS2)], colors{1}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '2hr +/- std');
plot(time2, mfitS2, 'Color', colors{1}, 'LineStyle', lineStyles{1}, 'DisplayName', '2hr @ 27C');
fill([time2, fliplr(time2)], [mfitT2 + stdT2, fliplr(mfitT2 - stdT2)], colors{1}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time2, mfitT2, 'Color', colors{1}, 'LineStyle', lineStyles{2}, 'DisplayName', '2hr @ 37C');

legend('Location', 'southwest');
title('2 Hr Sample Group Means', 'FontWeight','normal');
xlabel('Time (min)'); ylabel('Normalized Area');
xlim([0, 125]); xticks(linspace(0, 120, 9));
ylim([0.8, 1]); yticks(linspace(0.8, 1, 5));
pbaspect([1 1 1]); grid off; hold off;

% --- Standalone 4hr Figure ---
fig4 = figure('Units', 'inches', 'Position', [1, 1, height * 2, height]);
hold on;
fill([time4, fliplr(time4)], [mfitS4 + stdS4, fliplr(mfitS4 - stdS4)], colors{2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '4hr +/- std');
plot(time4, mfitS4, 'Color', colors{2}, 'LineStyle', lineStyles{1}, 'DisplayName', '4hr @ 27C');
fill([time4, fliplr(time4)], [mfitT4 + stdT4, fliplr(mfitT4 - stdT4)], colors{2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time4, mfitT4, 'Color', colors{2}, 'LineStyle', lineStyles{2}, 'DisplayName', '4hr @ 37C');

legend('Location', 'southwest');
title('4 Hr Sample Group Means', 'FontWeight','normal');
xlabel('Time (min)'); ylabel('Normalized Area');
xlim([0, 245]); xticks(linspace(0, 240, 17));
ylim([0.8, 1]); yticks(linspace(0.8, 1, 5));
pbaspect([2 1 1]); grid off; hold off;


% --- Standalone 6hr Figure ---
fig6 = figure('Units', 'inches', 'Position', [1, 1, height * 3, height]);
hold on;
fill([time6, fliplr(time6)], [mfitS6 + stdS6, fliplr(mfitS6 - stdS6)], colors{3}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '6hr +/- std');
plot(time6, mfitS6, 'Color', colors{3}, 'LineStyle', lineStyles{1}, 'DisplayName', '6hr @ 27C');
fill([time6, fliplr(time6)], [mfitT6 + stdT6, fliplr(mfitT6 - stdT6)], colors{3}, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(time6, mfitT6, 'Color', colors{3}, 'LineStyle', lineStyles{2}, 'DisplayName', '6hr @ 37C');

legend('Location', 'southwest');
title('6 Hr Sample Group Means', 'FontWeight','normal');
xlabel('Time (min)'); ylabel('Normalized Area');
xlim([0, 365]); xticks(linspace(0, 360, 25));
ylim([0.8, 1]); yticks(linspace(0.8, 1, 5));
pbaspect([3 1 1]); grid off; hold off;

% save all 

saveas(fig2, root + "/FIGURES/Normalized Area over Time (2Hr Sample Means).svg", 'svg');
saveas(fig4, root + "/FIGURES/Normalized Area over Time (4Hr Sample Means).svg", 'svg');
saveas(fig6, root + "/FIGURES/Normalized Area over Time (6Hr Sample Means).svg", 'svg');

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

%% PLOT ALL GROUPS ACROSS 6 TILES (DYNAMIC SWELLING)

% Bundle group data, labels, file suffixes, and time vectors into cell arrays for indexing
groupData   = {S2, S4, S6, T2, T4, T6};
groupMeans = {mS2; mS4; mS6; mT2; mT4; mT6};
colorshex2 = {'#5e6e3a', '#e99f78', '#f088f7', '#5e6e3a', '#e99f78', '#f088f7'};
groupLabels = {'2hr @ 27C', '4hr @ 27C', '6hr @ 27C', '2hr @ 37C', '4hr @ 37C', '6hr @ 37C'};
fileSuffix  = {'S2', 'S4', 'S6', 'T2', 'T4', 'T6'};
timeVectors = {time2, time4, time6, time2, time4, time6};

% Initialize the master canvas layout (2 rows x 3 columns)
masterFig = figure('Units', 'normalized', 'Position', [0.05, 0.05, 0.90, 0.85]);
layoutGrid = tiledlayout(masterFig, 2, 3, 'TileSpacing', 'loose', 'Padding', 'loose');

% Global structural annotations across layout borders
ylabel(layoutGrid, 'Normalized Area', 'FontSize', 13, 'FontWeight', 'bold');
xlabel(layoutGrid, 'Time (min)', 'FontSize', 13, 'FontWeight', 'bold');
title(layoutGrid, 'Exponential Fits for Normalized Area over Time (Individual Replicates)', 'FontSize', 15, 'FontWeight', 'bold');

% Iterate over the 6 unique experimental conditions
for g = 1:6
    % Extract localized configurations for this single tile index
    currentData  = groupData{g};
    currentMean = groupMeans{g};
    currentLabel = groupLabels{g};
    currentTime  = timeVectors{g};
    numPts       = length(currentTime);
    
    % Initialize storage arrays for metrics tracking
    fitG = zeros(rep, numPts);
    Rs   = zeros(rep, 1);
    fitGM = zeros(rep, numPts);
    RsM   = zeros(rep, 1);
    
    % Step into the target layout slot
    ax = nexttile(layoutGrid);
    hold(ax, 'on');
    
    % Process and display each discrete replicate sample
    for i = 1:rep
        % Calculate isolated fit metrics cleanly anchored to 1
        [fG, gof] = fit(currentTime', currentData(i, :)', ft, opts);
        Rs(i) = gof.rsquare;
        fitG(i, :) = fG(currentTime')';

        
        % Plot raw components and overlay fit line cleanly (suppressing redundant legends)
        plot(ax, currentTime, fitG(i, :), 'Color', colorshex2{g}, 'HandleVisibility', 'off');
        plot(ax, currentTime, currentData(i, :), 'Color', 'k', 'Marker', 'o', 'LineStyle', '--', 'HandleVisibility', 'off');

        
    end

    [fGM, gofM] = fit(currentTime', currentMean', ft, opts);
    RsM(g) = gofM.rsquare;
    fitGM(g, :) = fGM(currentTime')';

    plot(ax, currentTime, fitGM(g, :), 'Color', colorshex2{g}, 'LineWidth', 2.5, 'HandleVisibility', 'off');

    % =========================================================================
    % CALCULATE THE LOW GLOBAL R^2 AGAINST ALL RAW REPLICATE DATA POINTS
    % =========================================================================
    % 1. Flatten the raw data matrix row-by-row into a single column vector
    y_raw = currentData'; 
    y_raw = y_raw(:); 
    
    % 2. Replicate the time column vector to match the total pooled point count
    t_extended = repmat(currentTime', rep, 1);
    
    % 3. Project predictions from the mean model over the long time axis
    y_pred = fGM(t_extended);
    
    % 4. Dynamically strip out NaN values safely before manual R^2 calculations
    nanIdx = isnan(y_raw) | isnan(y_pred);
    y_raw(nanIdx) = [];
    y_pred(nanIdx) = [];
    
    % 5. Manual baseline calculation vs all combined variance noise
    ss_res = sum((y_raw - y_pred).^2);
    ss_tot = sum((y_raw - mean(y_raw)).^2);
    RsM(g) = 1 - (ss_res / ss_tot);

    
    % Generate dummy lines to isolate clean global legend references per tile
    plot(ax, NaN, NaN, 'Color', 'k', 'Marker', 'o', 'LineStyle', '--', 'DisplayName', currentLabel);

    plot(ax, NaN, NaN, 'Color', colorshex2{g}, 'DisplayName', 'individual samples');
    plot(ax, NaN, NaN, 'Color', colorshex2{g}, 'LineWidth', 2.5, 'DisplayName', 'group mean');
    
    % Customize axis boundaries, formatting rules, and labels for this specific group
    title(ax, currentLabel, 'FontWeight', 'bold', 'FontSize', 11);
    legend(ax, 'show', 'Location', 'northeast', 'FontSize', 9);
    
    % Enforce user constraints directly
    xlim(ax, [0, max(currentTime) + 5]);
    %ylim(ax, [0.75, 1.0]);
    
    grid(ax, 'off');
    hold(ax, 'off');

    r2string = sprintf('Samples Mean R^2: %.4f\n', mean(Rs));
    fprintf(r2string);
    r2Mstring = sprintf('Group Mean R^2: %.4f\n', mean(RsM(g)));
    fprintf(r2Mstring);
    
end

% Save out the composite dashboard visualization
%saveas(masterFig, root + "/FIGURES/Normalized_Area_All_Groups_Grid.svg", 'svg');
