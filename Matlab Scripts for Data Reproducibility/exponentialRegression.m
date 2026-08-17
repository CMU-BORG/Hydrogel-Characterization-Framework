%% initialize data arrays
clear all; close all; %clc

currentFolder = pwd;
root = currentFolder + "/Experimental Data/";

inventory_file = "Inventory Final.xlsx";
inventory = readtable(root + inventory_file,"ReadVariableNames",true);

inventoryE = inventory(:, [3, 4, 16]);

% Sort data by crosslinkTemp
sortedDataE = sortrows(inventoryE, 'CrosslinkTemp');

% Separate data by temperature
data27E = sortedDataE(sortedDataE.CrosslinkTemp == 27, :);
data37E = sortedDataE(sortedDataE.CrosslinkTemp == 37, :);

% Include entries where CrosslinkTime is 0 in both datasets
data27E = [data27E; inventoryE(inventoryE.CrosslinkTime == 0, :)];
data37E = [data37E; inventoryE(inventoryE.CrosslinkTime == 0, :)];

% Remove duplicate rows
data27E = unique(data27E, 'rows');
data37E = unique(data37E, 'rows');


% Combine data and create the X variable
combinedData = [data27E; data37E];
X = [zeros(height(data27E), 1); ones(height(data37E), 1)];
t = combinedData.CrosslinkTime;
E = combinedData.EffectiveElasticModulus;

%% DYNAMIC SWELLING EXPONENTIAL REGRESSION ANALYSIS

A_file = "Areas Final new.csv";
A = readmatrix(root+A_file);

Anorm = A(:, :)./A(:, 1) * 100;

Anorm(Anorm==0) = nan;

rows27 = find(inventory.CrosslinkTemp == 27)';
rows37 = find(inventory.CrosslinkTemp == 37)';

% Time intervals (15 min intervals from t=0 to t=360)
time_intervals = 0:15:360;
time_groups = [2,4,6];

% Storage for coefficients
coefficients_by_time = struct();

% Loop through each time group
for tg = time_groups

    fprintf("\n=== TIME GROUP: %dh ===\n", tg);

    % Rows belonging to this time group
    tg_rows_27 = rows27(inventory.CrosslinkTime(rows27) == tg);
    tg_rows_37 = rows37(inventory.CrosslinkTime(rows37) == tg);

    % Build dataset for this time group only
    data = [];

    % 27°C entries
    for row = tg_rows_27
        for t = time_intervals
            val = Anorm(row, t/15 + 1);
            if ~isnan(val)
                data = [data; t, 0, val];   % X = 0 for 27°C
            end
        end
    end

    % 37°C entries
    for row = tg_rows_37
        for t = time_intervals
            val = Anorm(row, t/15 + 1);
            if ~isnan(val)
                data = [data; t, 1, val];   % X = 1 for 37°C
            end
        end
    end

    % Convert to table
    data_table = array2table(data, 'VariableNames', {'t', 'X', 'Anorm'});
    
    % Define the fit equation
    fit_func = @(b, t, X) b(1) * exp(-t .* (b(2) + b(3) .* X)) + b(4) + b(5) .* X;
    
    % Extract t, X, and Anorm values for fitting
    t_values = data_table.t;
    X_values = data_table.X;
    Anorm_values = data_table.Anorm;
    
    % Initial guess for the coefficients
    initial_guess = [15, .1, .1, 80, 1];
    
    % Define lower bounds for the coefficients
    lower_bounds = [0, 0, 0, 0, 0];
    
    % Fit the nonlinear model with bounds
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    [b, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(@(b, tX) fit_func(b, tX(:,1), tX(:,2)), initial_guess, [t_values, X_values], Anorm_values, lower_bounds, [], options);
    
    % Store coefficients
    coefficients_by_time.(sprintf("t%dh", tg)) = b;

    % Display results
    fprintf("Fitted coefficients for %dh:\n", tg);
    disp(b);

    % R-squared
    SS_res = sum(residual.^2);
    SS_tot = sum((Anorm_values - mean(Anorm_values)).^2);
    R_squared = 1 - SS_res / SS_tot;

    fprintf("R-squared for %dh: %.4f\n", tg, R_squared);

    % Calculate statistics
    ci = nlparci(b, residual, 'jacobian', jacobian)
    se = sqrt(diag(inv(jacobian' * jacobian) * resnorm / (length(Anorm_values) - length(b))));
    t_stats = b ./ se;
    p_values = 2 * (1 - tcdf(abs(t_stats), length(Anorm_values) - length(b)));

end

% Final summary
disp('=======================================');
disp('SUMMARY OF FITTED COEFFICIENTS BY TIME');
disp(coefficients_by_time);



%% FUNCTIONS

function summary = nonlinearRegressionSummary(fit_func, t, X, y, initialGuess)
% nonlinearRegressionSummary - Fits a nonlinear model and returns a summary
%
% Inputs:
%   fit_func     - function handle: @(b, t, X) ...
%   t            - vector of independent variable
%   X            - vector of second predictor (e.g., group)
%   y            - response variable
%   initialGuess - initial guess for coefficients
%
% Output:
%   summary - struct containing coefficients, SE, tStat, pValue, CI, R2, R2_adj, MSE

    % Combine predictors
    xData = [t, X];

    % Fit model
    [b, R, J, CovB, mse] = nlinfit(xData, y, @(b, x) fit_func(b, x(:,1), x(:,2)), initialGuess);

    % Standard errors
    SE = sqrt(diag(CovB));

    % t-statistics
    tStat = b ./ SE;

    % Degrees of freedom
    df = length(y) - length(b);

    % p-values
    pValue = 2 * (1 - tcdf(abs(tStat), df));

    % Confidence intervals
    CI = nlparci(b, R, 'covar', CovB);

    % Fitted values
    yFit = fit_func(b, t, X);

    % R-squared
    SStot = sum((y - mean(y)).^2);
    SSres = sum((y - yFit).^2);
    R2 = 1 - SSres / SStot;
    R2_adj = 1 - (1 - R2) * (length(y) - 1) / df;

    % Package results
    summary = struct();
    summary.Coefficients = b;
    summary.SE = SE;
    summary.tStat = tStat;
    summary.pValue = pValue;
    summary.CI = CI;
    summary.R2 = R2;
    summary.R2_adj = R2_adj;
    summary.MSE = mse;

    % Display summary
    fprintf('\n%-10s %-12s %-12s %-12s %-12s\n', 'Coeff', 'Estimate', 'SE', 'tStat', 'pValue');
    for i = 1:length(b)
        fprintf('%-10s %12.4f %12.4f %12.4f %12.4g\n', ...
            ['b' num2str(i)], b(i), SE(i), tStat(i), pValue(i));
    end
    fprintf('\nR-squared: %.4f\n', R2);
    fprintf('Adjusted R-squared: %.4f\n', R2_adj);
    fprintf('MSE: %.4f\n', mse);
end
