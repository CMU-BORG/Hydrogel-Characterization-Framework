%% Fitting Exponential Models for Dynamic Area Data to Obtain Sample Parameters
%{

%}

clear all; close all; clc

%% Loading the Experimental data

% data location
% currentFolder = pwd;
root = "../Experimental Data/";

% loading the inventory table with all of the summary data
inventory_file = "Inventory Final.xlsx";
inventory = readtable(root + inventory_file,"ReadVariableNames",true);

% pulling out relevant columns of the inventory table (CrosslinkTime,
% CrosslinkTemp, Scale_pxpmm)
inventoryE = inventory(:, [3, 4, 12]);

% Loading the cross-sectional area data
A_file = "Areas Final new.csv";
A = readmatrix(root+A_file);    % [mm^2] CSA at each frame (15 min interval)
A(A==0) = nan;

rows27 = find(inventory.CrosslinkTemp == 27)';
rows37 = find(inventory.CrosslinkTemp == 37)';
rows = [rows27, rows37];

% Time intervals (15 min intervals from t=0 to t=360)
time_intervals = 0:15:360;

A_func = @(t,theta) (theta(1) + theta(2)) - theta(2)*exp(-t/theta(3));

% preallocating memory for the fitted exponential models (X, A0, deltaA,
% tau), where X=0 for low temp (27C) and X=1 for high temp (37C)
params = zeros(length(rows27)+length(rows27),4);

f27 = figure('Position',[100,100,600,300],'Color','w'); hold all
f37 = figure('Position',[700,100,600,300],'Color','w'); hold all

for i = 1:length(rows)
    r = rows(i);
    temp_i = inventoryE{r,"CrosslinkTemp"}; % [C]

    A_i = A(r,:); % getting the area data for given sample
    t_i = time_intervals;   % time vector
    t_i(isnan(A_i)) = [];   % remove the nan values to prevent interferring with fitting
    A_i(isnan(A_i)) = [];

    err = @(params) sum( (A_i - A_func(t_i,params)).^2 );
    params_i = fminsearch(err,[A_i(1),A_i(end)-A_i(1),15]);
    
    % compiling results
    X_i = (temp_i - 27)/10; % difference in temperature from the control group (27 C)
    params(i,:) = [X_i, params_i];
    
    if X_i==0
        figure(f27); hold all
        title("T=27^{\circ}C")
        col = [1,0.25,0.25]*(0.5 + 0.5*rand(1));
    else
        figure(f37); hold all
        title("T=37^{\circ}C")
        col = [0.25,0.25,1]*(0.5 + 0.5*rand(1));
    end
    plot(t_i,A_i,'ok','MarkerFaceColor',col)
    plot(t_i,A_func(t_i,params_i),'-','color',col,'LineWidth',2)
    xlabel("Time [min]")
    ylabel("Cross-sectional Area [mm^2]")
    xlim([0,60*6])
    ylim([100,275])
end

A0_27 = params(params(:,1)==0,2);
dA_27 = params(params(:,1)==0,3);
tau_27 = params(params(:,1)==0,4);

A0_37 = params(params(:,1)==1,2);
dA_37 = params(params(:,1)==1,3);
tau_37 = params(params(:,1)==1,4);

fprintf("~~~~ T=27 C ~~~~\n");
fprintf("\tA_0 = %.3f +/- %.3f\n",mean(A0_27),std(A0_27))
fprintf("\tdA = %.3f +/- %.3f\n",mean(dA_27),std(dA_27))
fprintf("\ttau = %.3f +/- %.3f\n",mean(tau_27),std(tau_27))

fprintf("~~~~ T=37 C ~~~~\n");
fprintf("\tA_0 = %.3f +/- %.3f\n",mean(A0_37),std(A0_37))
fprintf("\tdA = %.3f +/- %.3f\n",mean(dA_37),std(dA_37))
fprintf("\ttau = %.3f +/- %.3f\n",mean(tau_37),std(tau_37))

save("AreaModelParameters.mat","params",'-mat')