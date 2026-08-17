%% calculating monomer ratios

clear all; close all; %clc

currentFolder = pwd;
root = currentFolder + "/Experimental Data/";

filename = "mldcon.csv";
sample1 = readtable(root + "/Alginate1-10/10/pdata/1/" + filename, "ReadVariableNames",true);
sample2 = readtable(root + "/Alginate2-60Deg/10/pdata/1/" + filename, "ReadVariableNames",true);
sample3 = readtable(root + "/Alginate3-60Deg/10/pdata/1/" + filename, "ReadVariableNames",true);

% peak assignments for each sample
peaks1 = [5.361523061, 5.076114542, 5.03800706, 5.002526049, 4.96568006, 4.765764841];
peaks2 = [5.338987188, 5.055231478, 5.016641198, 4.979523036, 4.943317431, 4.741294078];
peaks3 = [5.338626538, 5.055360183, 5.016879853, 4.980163016, 4.941841015, 4.740802035];

integrals = zeros(3, length(peaks1));
integrals2 = size(length(peaks2));
integrals3 = size(length(peaks3));

for i = 1:length(peaks1)

integrals(1, i) = sample1.Area_abs_Hz_(round(sample1.Frequency_ppm_, 4) == round(peaks1(i), 4));
integrals(2, i) = sample2.Area_abs_Hz_(round(sample2.Frequency_ppm_, 4) == round(peaks2(i), 4));
integrals(3, i) = sample3.Area_abs_Hz_(round(sample3.Frequency_ppm_, 4) == round(peaks3(i), 4));

end

AA = integrals(:, 1);
AB = integrals(:, 2) + integrals(:, 3) + integrals(:, 4) + integrals(:, 5);
AC = integrals(:, 6);

Gfrac = AA ./ (AB + AC)

MGRatio = (1 - Gfrac)./Gfrac

Gdoublet = AC ./ (AB + AC)
Mdoublet = abs((1 - Gfrac) - (Gfrac - Gdoublet))

MGDoubletRatio = Mdoublet ./ Gdoublet


% % Perform one-way ANOVA
% [p, tbl, stats] = anova1(MGRatio);
% 
% % Display the p-value
% disp(['p-value from ANOVA: ', num2str(p)]);
% 
% % Perform pairwise t-tests
% [h12, p12] = ttest2(MGRatio(1), MGRatio(2));
% [h13, p13] = ttest2(MGRatio(1), MGRatio(3));
% [h23, p23] = ttest2(MGRatio(2), MGRatio(3));
% 
% % Display the p-values
% disp(['p-value for value1 vs value2: ', num2str(p12)]);
% disp(['p-value for value1 vs value3: ', num2str(p13)]);
% disp(['p-value for value2 vs value3: ', num2str(p23)]);
