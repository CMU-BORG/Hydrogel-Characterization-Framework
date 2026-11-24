%%  ABC-SMC Sampler for Hydrogel Stiffness Progression Model
%{
    This script will perform the Approximate Bayesian Computation (ABC)
    Sequential Monte Carlo (SMC) calibration of the statistical parameters
    governing the hydrogel stiffness progression model. For a given sample
    of the hydrogel (i), it is assumed that its stiffness increases at a
    decreasing rate from an initial condition (E0) to a steady state 
    stiffness (E0 + dE) with a particular time constant (tau):
        
            E_i(t) = (E0_i + dE_i) - E0_i*exp(-t / tau)

    We assume that the values of E0, dE, and tau are not constant but are
    different for each sample due to differences in composition, exact
    fabrication condition, etc. Therefore, we construct a probabilistic
    model for the distribution of E0, dE, and tau such that we can sample
    from this model to generate realistic models of hydrogel samples.

    For this statistical model, we assume that initial stiffness and final
    stiffness are drawn from a joint normal distribution (capturing the
    fact that stiffness is related to monomer composition, so compositions
    that start stiffer are likely to end up stiffer). We chose normal or
    uniform distributions to keep priors uninformed and to prevent 
    potential skewness of the distribution. We assume the time constant is 
    independent of the stiffness values and is modeled as a Gamma random 
    variable. 

    This combined statistical model has 7 parameters to estimate: mu and
    sigma for each of the normal marginal distributions (4), the
    correlation coefficient of the normal distribution (1), and alpha
    and beta of the gamma distribution (2). These parameters will be
    calibrated using the ABC-SMC algorithm.

    This same procedure will work for any statistical model you want to
    try. The specifics (what parameters you're dealing with, the priors,
    the specific stats model, and the simulated experiments) will change,
    but the overall flow will be the same.
%}

clear all; close all; clc;

% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool("local", 8);  % Adjust the number of workers if needed
end

timestamp = datestr(now, 'yyyymmdd_HHMM');

%% Loading experimental data and calculating the summary metrics for each time point in the experiment (FILL IN)

% load the data that you have from your experiments
root = "/placeholder"
inventory_file = "placeholder.xlsx";
inventory = readtable(root + inventory_file,"placeholder",true);

% extract relevant sample metrics
inventory = inventory(:, [columns]);

% Sort data by desired metric
sortedData = sortrows(inventory, 'metric');
dataSubset = sortedData(sortedData.metric == X, :);

% define time groups
time1 = dataSubset.metric(dataSubset.time == 1)';
time2 = dataSubset.metric(dataSubset.time == 2)';
time3 = dataSubset.metric(dataSubset.time == 3)';
time4 = dataSubset.metric(dataSubset.time == 4)';

timedata = [time1; time2; time3; time4];
time_points = [1, 2, 3, 4];  

% For each time point, you will calculate the summary metrics that will be
% compared to your simulations during each step of the tuning. You will
% want a metric to capture the location (central tendency) of the data
% (mean, median, etc) and to capture the spread of the data (std,
% quartiles, etc). There isn't necessarily a standard for this, but if your
% data is non-normal, then something like median and the IQR would make
% sense.

N_timepoints = 4; 

sum_met_central = zeros(N_timepoints,1); % summary metric of central tendency for each time point (mean, median, etc)
sum_met_spread = zeros(N_timepoints,1);     % summary metric of spread for each time point (std, quartiles, etc)
samp_size = zeros(N_timepoints,1);          % number of samples at each time point

f_sum_met_central = @(x) median(x);    % add whatever central tendency metric you want to use. doing it programmatrically 
f_sum_met_spread = @(x) iqr(x);    % like this will make it easier to change things later and you will use them in the simulations later

for index=1:N_timepoints
    x_i = timedata(index, :); % pull out the data for just this time point
    sum_met_central(index) = f_sum_met_central(x_i); 
    sum_met_spread(index) = f_sum_met_spread(x_i);
    samp_size(index) = length(x_i);
end
N_full = sum(samp_size);    % total number of samples to generate per experiment


%% Setting up the prior distributions on the statistical parameters (FILL IN)

%{
NOTES ON PARAMETER PRIORS: 

Because there isn't much existing information about hydrogels and the 
timecourse of their stiffening, it makes sense to use fairly uninformative 
priors. That way most of the information we get from the final parameters 
will be controlled by the data, not by the priors.

Because we assume the initial stiffness (E0) and stiffness chnge (dE) to be
within the same order of magnitude based on literature, we will use the 
same priors for the parameters of their distributions. 

The distribution is parameterized by mu (median)and sigma (standard 
deviation). mu can be anyvalue, but sigma must be positive. We can use a 
normal distribution for mu and a uniform distribution for sigma.
%}

mu_E0_mean = ; % [metric unit] (order of magnitude of metric)
mu_E0_std = ;  % [metric unit] (assume 95% of prior distribution is between _ and _)
mu_E0_prior = makedist("Normal","mu",mu_E0_mean, "sigma",mu_E0_std);


mu_dE_mean = ; % [metric unit] (order of magnitude of metric)
mu_dE_std = ;  % [metric unit] (assume 95% of prior distribution is between _ and _)
mu_dE_prior = makedist("Normal","mu",mu_dE_mean, "sigma",mu_dE_std);

sig_E_lower = ;      % [metric unit] lowest possible value of sigma (E within +/- 0.3% of the mean) 
sig_E_upper = ;    % [metric unit] highest possible values (needs to be tuned)
sig_E_prior = makedist("Uniform","lower",sig_E_lower,"upper",sig_E_upper);

% CORRELATION COEFFICIENT PRIOR: The correlation coefficient of the
% multivariate normal is constrained to  [-1,1]. Because we don't have
% information about what it should be we can just use a uniform
% distribution prior.
rho_E_prior = makedist("Uniform", "lower", -1, "upper", 1);

% GAMMA PARAMETER PRIORS: The gamma distribution is parameterized by alpha
% and beta, but it can reparameterized by the mean and standard deviation 
% (mean = alpha*beta, variance = beta*mean). Both of these parameters have
% to be positive, so it would make sense to put gamma distribution priors
% on both. Alternatively you could use a uniform distribution on the either
% or both

% Prior for mean of tau 
mean_tau_prior = makedist("Gamma", "a", , "b", );  % [time unit]
sig_tau_prior = makedist("Gamma", "a", , "b", );  % [time unit]

%% Sequential Monte Carlo Approximate Bayesian Computation (FILL IN)
 
N_particles = 1000;     % number of particles in each approximate distribution

% Specify the width of the random walk that the particles will
% take. This might take a little bit of tuning to make sure the sampler
% converges well. 
meanE_proposal_width = ;   % [metric unit]
stdE_proposal_width = ;    % [metric unit]
rho_proposal_width = ;   % [  ]
meanTau_proposal_width = ; % [time unit]
sigTau_proposal_width = ;  % [time unit]

T = 20;  % number of iteration. this can also be played with until the solver converges well

% During the SMC, decrease the acceptable error between the
% simulations and the data. This version uses an exponential interpolation
% to reduce the step size as the acceptance threshold decreases to help
% with efficient convergence.

% epsilon as percent diff 
epsilon_sum_met_central = exp(linspace(log(100), log(5), T));  % units: kPa 
epsilon_sum_met_spread = exp(linspace(log(100), log(10), T));  % units: kPa

% variables to store {weights, mean_E0, sig_E0, mean_Ef, sig_Ef, rho_samps, mean_tau, sig_tau, samples_parameters}
intermediate_distributions = cell(T,9);     % storing the intermediate distributions for later visualization

% preallocating memory -- these will be reused for each population
mean_E0_samps = zeros(N_particles,1);
sig_E0_samps = zeros(N_particles,1);
mean_dE_samps = zeros(N_particles,1);
sig_dE_samps = zeros(N_particles,1);
rho_samps = zeros(N_particles,1);
mean_tau_samps = zeros(N_particles,1);
sig_tau_samps = zeros(N_particles,1);


% weight vector that is used to weight the sampling of the intermediate
% distributions
weights = zeros(N_particles,1);
indexes = 1:N_particles;

%% Running the Sampler
for kk = 1:T
    % for each intermediate distribution

    if kk>1
        % pulling previous intermediate distribution. Except in the first
        % round, you will sample from these distributions
        weights_int     = intermediate_distributions{kk-1,1};
        mean_E0_int     = intermediate_distributions{kk-1,2};
        sig_E0_int      = intermediate_distributions{kk-1,3};
        mean_dE_int     = intermediate_distributions{kk-1,4};
        sig_dE_int      = intermediate_distributions{kk-1,5};
        rho_int         = intermediate_distributions{kk-1,6};
        mean_tau_int    = intermediate_distributions{kk-1,7};
        sig_tau_int     = intermediate_distributions{kk-1,8};
    end

    samples_parameters = 0;
    std_fails = 0;
    mean_fails = 0;

    accepted_particles = 0;

    parfor index = 1:N_particles %parfor
        % loop through all of your particles
        accept = 0;

        while ~accept
            % keep trying until you get a parameter that satisfies the current error bound
            samples_parameters = samples_parameters + 1; % keep track of total number of sims needed

            % sample parameters from intermediate distribution with weights (if
            % kk==1, sample from prior)
            if kk==1
                % sample values from your priors
                mean_E0_prop = mu_E0_prior.random(1); 
                mean_dE_prop = mu_dE_prior.random(1); 

                sig_E0_prop = sig_E_prior.random(1);
                sig_dE_prop = sig_E_prior.random(1);

                rho_prop = rho_E_prior.random(1);

                mean_tau_prop = mean_tau_prior.random(1);
                sig_tau_prop = sig_tau_prior.random(1);

            else
                % randomly sample a particle index weighted by the previous
                % interations weights. Here, you sample an index so you
                % make sure that you pull all of the parameters from a
                % single particle so that you properly maintain the
                % correlations between steps.
                prop_index = randsample(indexes,1,true,weights_int);

                % get a value from the previous distribution
                mean_E0_cur = mean_E0_int(prop_index);
                sig_E0_cur = sig_E0_int(prop_index);
                mean_dE_cur = mean_dE_int(prop_index);
                sig_dE_cur = sig_dE_int(prop_index);
                rho_cur = rho_int(prop_index);

                mean_tau_cur = mean_tau_int(prop_index);
                sig_tau_cur = sig_tau_int(prop_index);
                
                % ~~~ perturb the resampled parameters. As a first pass, I 
                % would just perturb them with a random walk (add on Gaussian
                % noise with the widths laid out above. You will need to make
                % sure that the values stay within the bounds. ~~~

                mean_E0_prop = mean_E0_cur + normrnd(0,meanE_proposal_width);
                mean_dE_prop = mean_dE_cur + normrnd(0, meanE_proposal_width);

                % check if bounded variables are within bounds. If not,
                % resample until they are

                sig_E0_prop = sig_E0_cur + normrnd(0, stdE_proposal_width);
                while sig_E0_prop <= 0
                    sig_E0_prop = sig_E0_cur + normrnd(0, stdE_proposal_width);
                end                
                
                sig_dE_prop = sig_dE_cur + normrnd(0, stdE_proposal_width);
                while sig_dE_prop <= 0
                    sig_dE_prop = sig_dE_cur + normrnd(0, stdE_proposal_width);
                end

                rho_prop = rho_cur + normrnd(0, rho_proposal_width);
                while rho_prop < -1 || rho_prop > 1
                    rho_prop = rho_cur + normrnd(0, rho_proposal_width);
                end

                mean_tau_prop = mean_tau_cur + normrnd(0, meanTau_proposal_width);
                while mean_tau_prop <= 0
                    mean_tau_prop = mean_tau_cur + normrnd(0, meanTau_proposal_width);
                end

                sig_tau_prop = sig_tau_cur + normrnd(0, sigTau_proposal_width);
                while sig_tau_prop <= 0
                    sig_tau_prop = sig_tau_cur + normrnd(0, sigTau_proposal_width);
                end

            end

            % RUNNING THE SIMULATIONS: To conduct your simulated
            % experiments, you will first build you statistical model for
            % the proposed parameters. Then for each time point, you will
            % sample n_i values of E0, dE, and tau, and calculate the
            % stiffness of that sample based on the stiffening model. 
            % Finally, for each time point, you will use the summary metric 
            % functions you defined above to get the summary metrics from 
            % your simulation.
            
            sum_met_central_sim = zeros(N_timepoints,1); 
            sum_met_spread_sim = zeros(N_timepoints,1);  

            % build your statistical model. For the multivariate lognormal,
            % you will need to construct the covariance matrix from the
            % proposed sig_E0, sig_dE, and rho. 
            
            mu_E = [mean_E0_prop, mean_dE_prop];
            cov_E = [sig_E0_prop^2, rho_prop*sig_E0_prop*sig_dE_prop; rho_prop*sig_E0_prop*sig_dE_prop, sig_dE_prop^2];
            E_dist = @(n) mvnrnd(mu_E, cov_E,n) ; % samples n pairs of values from a multivariate norm

            gam_alpha = (mean_tau_prop / sig_tau_prop)^2;
            gam_beta = sig_tau_prop^2 / mean_tau_prop;
            tau_dist = @(n) gamrnd(gam_alpha,gam_beta,n, 1);
            
            for counter = 1:N_timepoints
                
                n = samp_size(counter);

                % get n samples of [E0,dE] and tau
                E0_dE = E_dist(n);         % n x 2 matrix: [E0, dE]
                tau = tau_dist(n);         % n x 1 vector

                t = time_points(counter);        % current time point

                % for each set of parameters calculate the stiffness based
                % on the stiffness model

                E = (E0_dE(:,1) + E0_dE(:,2)) - E0_dE(:,2) .* exp(-t ./ tau);
                
                % calculate the summary metrics
                sum_met_central_sim(counter) = f_sum_met_central(E); 
                sum_met_spread_sim(counter) = f_sum_met_spread(E);
            end


            % determining if you simulated data are close enough to the
            % experimental data.

            % epsilon as a percent diff
            err_central = 100 * sqrt( mean(( (sum_met_central - sum_met_central_sim) ./ sum_met_central ).^2));
            err_spread  = 100 * sqrt( mean(( (sum_met_spread - sum_met_spread_sim) ./ sum_met_spread ).^2));

            % check if the data is close enough in both central tendency
            % and spread
            accept_central =  (err_central < epsilon_sum_met_central(kk)); 
            accept_spread = (err_spread < epsilon_sum_met_spread(kk)); 

            % check distance to particle acceptance
            % print diagnostics
            disp("Generation " + kk + ", Particle " + index);
            disp("  err_central = " + string(err_central) + " vs epsilon = " + string(epsilon_sum_met_central(kk)));
            disp("  err_spread  = " + string(err_spread)  + " vs epsilon = " + string(epsilon_sum_met_spread(kk)));

            if ~accept_central
                disp("Rejected due to central error");
            end
            if ~accept_spread
                disp("Rejected due to spread error");
            end

            % determining if to break out of the while loop
            accept = accept_central & accept_spread;

            if accept
                accepted_particles = accepted_particles + 1;
            end
        end


        % you have now accepted this particle and it will be passed onto
        % the next round
        mean_E0_samps(index) = mean_E0_prop;
        sig_E0_samps(index) = sig_E0_prop;
        mean_dE_samps(index) = mean_dE_prop;
        sig_dE_samps(index) = sig_dE_prop;
        rho_samps(index) = rho_prop;
        mean_tau_samps(index) = mean_tau_prop;
        sig_tau_samps(index) = sig_tau_prop;


        % calculating new weight for particle. This weighting function is
        % found in the Toni et al. 2009 paper (in the algorithm, step S2.2)
        if kk==1
            % for the first distribution, you will just give all values
            % equal weight
            weights(index) = 1;
        else
            % setting up perturbation kernels on the new parameters

            % based on the random walk you conducted on the particles, you
            % will need to calculate the probability of going from the new
            % accepted particle to any of the values in the last
            % intermediate distribution (basically how likely was it that
            % you got to this value given all of the values in the previous
            % distribution). For the Gaussian random walk, this will come
            % from the pdf of a normal centered at the accepted particle
            % with std equal to the proposal width. I would use makedist to
            % create functions for each of these random walk distributions
            % because you will need them again.
            mean_E0_kernel = makedist("Normal","mu",mean_E0_prop,"sigma",meanE_proposal_width);
            sig_E0_kernel  = makedist("Normal", "mu", sig_E0_prop,  "sigma", stdE_proposal_width);
            mean_dE_kernel = makedist("Normal", "mu", mean_dE_prop, "sigma", meanE_proposal_width);
            sig_dE_kernel  = makedist("Normal", "mu", sig_dE_prop,  "sigma", stdE_proposal_width);
            rho_kernel     = makedist("Normal", "mu", rho_prop,     "sigma", rho_proposal_width);
            mean_tau_kernel = makedist("Normal", "mu", mean_tau_prop, "sigma", meanTau_proposal_width);
            sig_tau_kernel  = makedist("Normal", "mu", sig_tau_prop,  "sigma", sigTau_proposal_width);
            
            % getting probabilities of transitions for each parameter
            p_mean_E0 = mean_E0_kernel.pdf(mean_E0_int);
            p_sig_E0  = sig_E0_kernel.pdf(sig_E0_int);
            p_mean_dE = mean_dE_kernel.pdf(mean_dE_int);
            p_sig_dE  = sig_dE_kernel.pdf(sig_dE_int);
            p_rho     = rho_kernel.pdf(rho_int);
            p_mean_tau = mean_tau_kernel.pdf(mean_tau_int);
            p_sig_tau  = sig_tau_kernel.pdf(sig_tau_int);   

            % Next, get the total probability of the transition. Because 
            % all of the values were perturbed independently, this is just
            % the product of all the individual transition probabilities
            p_params = p_mean_E0 .* p_sig_E0 .* p_mean_dE .* p_sig_dE .* p_rho .* p_mean_tau .* p_sig_tau;

            % Then, you will calculate the prior probability of the accepted
            % particle. Again because parameters all have indepedent
            % priors, this is just the product of the prior probability of
            % all of the individual accepted parameters.
            p_prior = mu_E0_prior.pdf(mean_E0_prop) * ...
                      sig_E_prior.pdf(sig_E0_prop) * ...
                      mu_dE_prior.pdf(mean_dE_prop) * ...
                      sig_E_prior.pdf(sig_dE_prop) * ...
                      rho_E_prior.pdf(rho_prop) * ...
                      mean_tau_prior.pdf(mean_tau_prop) * ...
                      sig_tau_prior.pdf(sig_tau_prop);

            % calculating the updated weight
            weights(index) = p_prior / sum( weights_int .* p_params );
                        
        end
    end

        disp("Accepted " + string(accepted_particles) + " particles in generation " + string(kk));

        % normalizing new weights for the particles
        weights = weights ./ sum(weights);

        % save the intermediate distribution for post analysis
        intermediate_distributions(kk,:) = {weights, mean_E0_samps, sig_E0_samps, ...
        mean_dE_samps, sig_dE_samps, rho_samps, mean_tau_samps, sig_tau_samps, samples_parameters};

        save(root + "placeholder/ABC_SMC_intermediate_" + timestamp + ".mat", "intermediate_distributions");

%% === Plot mean and spread of accepted particles (no new variables) ===
        
    % plot simulated curves based on particles for every generation 
    plot_median = zeros(length(intermediate_distributions{kk,2}), 4);
    plot_iqr = zeros(length(intermediate_distributions{kk,2}), 4);
    
    % randomly sample 100 as before
    % for each particle, build statistical model
    % generate syntetic alginate data (add more time points)
    for i = 1:length(intermediate_distributions{kk,2})
        mu_E = [intermediate_distributions{kk,2}(i), intermediate_distributions{kk,4}(i)];
        cov_E = [intermediate_distributions{kk,3}(i)^2, ...
            intermediate_distributions{kk,6}(i)*intermediate_distributions{kk,3}(i)*intermediate_distributions{kk,5}(i); ...
            intermediate_distributions{kk,6}(i)*intermediate_distributions{kk,3}(i)*intermediate_distributions{kk,5}(i), ...
            intermediate_distributions{kk,5}(i)^2];
        E = mvnrnd(mu_E, cov_E, 100);
        tau = gamrnd((intermediate_distributions{kk,7}(i)/intermediate_distributions{kk,8}(i))^2, ...
            intermediate_distributions{kk,8}(i)^2/intermediate_distributions{kk,7}(i), 100, 1);
        E_curves = (E(:,1) + E(:,2)) - E(:,2) .* exp(-time_points ./ tau);
    
        plot_median(i, :) = median(E_curves);
        plot_iqr(i, :) = iqr(E_curves);
    end
    
    
    simulatedCurves = figure;
    
    % === Subplot 1: Median stiffness ===
    subplot(2,1,1); hold on;
    plot(time_points, sum_met_central, 'k--', 'LineWidth', 2);  % experimental median
    for i = 1:length(plot_median)
        plot(time_points, plot_median(i, :), 'Color', rand(1,3));
    end
    ylim([50 200]);  % fixed y-axis range
    title("Gen " + kk + ": Median Stiffness");
    xlabel("Time (min)"); ylabel("Median (kPa)");
    
    
    % === Subplot 2: IQR of stiffness ===
    subplot(2,1,2); hold on;
    plot(time_points, sum_met_spread, 'k--', 'LineWidth', 2);  % experimental median
    for i = 1:length(plot_iqr)
        plot(time_points, plot_iqr(i, :), 'Color', rand(1,3));
    end
    ylim([0 120]);  % fixed y-axis range
    title("Gen " + kk + ": IQR of Stiffness");
    xlabel("Time (min)"); ylabel("IQR (kPa)");
    
    saveas(simulatedCurves, root + "placeholder/simulatedCurves_" + timestamp + "_" + kk + ".svg", 'svg');

end

%% Final resampling of parameter distribution
% One last step that you can do is do one final
% weighted resampling of your parameter distribution. What this will do is
% reweight all of the samples to be the same, so you won't have to do any
% weighted calculations, and you can just directly calculate things like
% mean, median, etc from the resampled particles.

% last succesful generation
T = 20;
N_particles = 1000;

% Final weighted resampling of particles from last generation
final_weights     = intermediate_distributions{T,1};
final_mean_E0     = intermediate_distributions{T,2};
final_sig_E0      = intermediate_distributions{T,3};
final_mean_dE     = intermediate_distributions{T,4};
final_sig_dE      = intermediate_distributions{T,5};
final_rho         = intermediate_distributions{T,6};
final_mean_tau    = intermediate_distributions{T,7};
final_sig_tau     = intermediate_distributions{T,8};

% Resample indices based on final weights
resample_indices = randsample(1:N_particles, N_particles, true, final_weights);

% Create unweighted posterior samples
resampled_mean_E0   = final_mean_E0(resample_indices);
resampled_sig_E0    = final_sig_E0(resample_indices);
resampled_mean_dE   = final_mean_dE(resample_indices);
resampled_sig_dE    = final_sig_dE(resample_indices);
resampled_rho       = final_rho(resample_indices);
resampled_mean_tau  = final_mean_tau(resample_indices);
resampled_sig_tau   = final_sig_tau(resample_indices);

final_posterior = struct( ...
    "mean_E0", resampled_mean_E0, ...
    "sig_E0", resampled_sig_E0, ...
    "mean_dE", resampled_mean_dE, ...
    "sig_dE", resampled_sig_dE, ...
    "rho", resampled_rho, ...
    "mean_tau", resampled_mean_tau, ...
    "sig_tau", resampled_sig_tau ...
);

%save(root + "placeholder/ABC_SMC_final_" + timestamp + ".mat", "final_posterior");

summary_stats = struct( ...
    "median_mean_E0", median(resampled_mean_E0), ...
    "std_rho", std(resampled_rho), ...
    "tau_credible_interval", quantile(resampled_mean_tau, [0.025, 0.975]) ...
);

%save(root + "placeholder/ABC_SMC_summary_" + timestamp + ".mat", "summary_stats");



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Congratulations! At this point you have, successfully used the ABC-SMC
% algorithm to get a distribution parameters for your statistical model.
% Next, you can find you mean parameters (as the weighted mean of your
% particles from the final stage). And using those mean parameters, you can
% plot out what the mean curve and the std envelope for your model to see
% how it compares to your data (basically rerun your simulated experiments 
% from earlier but using the mean or median parameters). You can also look 
% at the uncertainty in your parameters by plotting the posterior joint 
% distribution of your parameters or calculate whatever statistics you 
% want for you parameters (mean, std, credible intervals, etc).
%
% Also definitely want to save all of the intermediate distributions and
% the final distribution because this may take a while to rerun.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%% === Posterior predictive visualization ===

% visualize posterior parameters in each generation
for i = 2:8
    figure(); hold all
    for j=1:T
        color = 1 - (j-1)/40 * [1,1,1];
        part = intermediate_distributions{j,i};
        weights = intermediate_distributions{j,1};
        part_weighted = randsample(part,1000,true,weights);
        histogram(part_weighted, 'FaceColor',color)

   end
end

%% === Population predictive visualization ===
time_vector = linspace(0, 360, 360);

% Posterior mean fit curve (get theta)
% sample statistical model
E0_mean   = mean(final_posterior.mean_E0);
sig_E0_mean = mean(final_posterior.sig_E0);
dE_mean   = mean(final_posterior.mean_dE);
sig_dE_mean = mean(final_posterior.sig_dE);

rho_mean = mean(final_posterior.rho);

tau_mean  = mean(final_posterior.mean_tau);
sig_tau_mean = mean(final_posterior.sig_tau);

% CI for posteriors
E0_CI   = prctile(final_posterior.mean_E0, [2.5, 97.5], 1);
sig_E0_CI = prctile(final_posterior.sig_E0, [2.5, 97.5], 1);
dE_CI   = prctile(final_posterior.mean_dE, [2.5, 97.5], 1);
sig_dE_CI = prctile(final_posterior.sig_dE, [2.5, 97.5], 1);

rho_CI = prctile(final_posterior.rho, [2.5, 97.5], 1);

tau_CI  = prctile(final_posterior.mean_tau, [2.5, 97.5], 1);
sig_tau_CI = prctile(final_posterior.sig_tau, [2.5, 97.5], 1);

%% generate matrix
mu_Epop = [E0_mean, dE_mean];
cov_Epop = [sig_E0_mean^2, rho_mean*sig_E0_mean*sig_dE_mean; ...
    rho_mean*sig_E0_mean*sig_dE_mean, sig_dE_mean^2];

% estimate how alginate parameters are distributed in population
Epop = mvnrnd(mu_Epop, cov_Epop, 100);
taupop = gamrnd((tau_mean/sig_tau_mean)^2, ...
    sig_tau_mean^2/tau_mean, 100, 1);

population = (Epop(:,1) + Epop(:,2)) - Epop(:,2) .* exp(-time_vector ./ taupop); % this is a matrix, rows are prob different versions in the popoulation, column prob time points
median_pop = median(population, 1);
% (median across a single time point)

%% calculate intervals and CI bounds (EDIT FOR DESIRED STATS METRICS)

% take 95% credible interval of generated alginate population
ci95_pop = prctile(population, [2.5, 97.5], 1);
ci50_pop = prctile(population, [25, 75], 1);
fit_lower_CI = ci95_pop(1, :);
fit_upper_CI = ci95_pop(2, :);

% take std of generated alginate population
std_pop = std(population, 1);

% take iqr of generated alginate population
iqr_pop = iqr(population, 1);

%% plot raw data, posterior mean fit curve, CI interval (EDIT FOR DATA SIZES AND DESIRED METRICS)

variabilityModel = figure;
hold on;

% +/- std shading
fill([time_vector fliplr(time_vector)], [fit_lower_std fliplr(fit_upper_std)], ...[0.8 0.8 1], 'EdgeColor','none','FaceAlpha',0.4); % std band

% 95% CI shading
fill([time_vector fliplr(time_vector)], [fit_lower_CI fliplr(fit_upper_CI)], ...
     [0.8 0.8 1], 'EdgeColor','none','FaceAlpha',0.2); % 95% band

% plot all population curves in light grey
plot(time_vector, population(i, :), '-', 'Color', [0.5, 0.5, 0.5, 0.3]);
for i = 2:100
    plot(time_vector, population(i, :), '-', 'Color', [0.5, 0.5, 0.5, 0.3], 'HandleVisibility','off');
end


% plot iqr brackets
errorbar(time_vector(1), median_pop(1), iqr_pop(1)/2, 'Color', 'b', 'LineWidth', 1);
errorbar(time_vector(120), median_pop(120), iqr_pop(120)/2, 'Color', 'b', 'LineWidth', 1, 'HandleVisibility','off');
errorbar(time_vector(240), median_pop(240), iqr_pop(240)/2, 'Color', 'b', 'LineWidth', 1, 'HandleVisibility','off');
errorbar(time_vector(360), median_pop(360), iqr_pop(360)/2, 'Color', 'b', 'LineWidth', 1, 'HandleVisibility','off');

% raw data points
scatter(dataSubset.CrosslinkTime, dataSubset.EffectiveElasticModulus, 25, "filled", 'MarkerFaceColor',"#4DBEEE");


% plot median curve
plot(time_vector, median_pop, 'k-', 'LineWidth', 2);


ylim([50 280]);
title("placeholder");
xlabel("Time ()"); ylabel("placeholder ()");
legend('95% CI', 'Simulated Particles', 'IQR', 'Raw Data Points', 'Median Posterior Fit', 'Location', 'southeast');

timestamp = datestr(now, 'yyyymmdd_HHMM');
saveas(variabilityModel, root + "placeholder/variabilityModel_" + timestamp + ".svg", 'svg');


%% ---- explicit posterior parameter output ----

Ef_samps = final_posterior.mean_E0 + final_posterior.mean_dE; % E0 + dE

fprintf('Posterior mu_E0: mean = %.2f, 95%% CI = [%.2f, %.2f]\n', ...
    mean(final_posterior.mean_E0), prctile(final_mean_E0,2.5), prctile(final_mean_E0,97.5));

fprintf('Posterior sigma_E0: mean = %.2f, 95%% CI = [%.2f, %.2f]\n', ...
    mean(final_posterior.sig_E0), prctile(final_sig_E0,2.5), prctile(final_sig_E0,97.5));

fprintf('Posterior mu_dE: mean = %.2f, 95%% CI = [%.2f, %.2f]\n', ...
    mean(final_posterior.mean_dE), ...
    prctile(final_mean_dE,2.5), ...
    prctile(final_mean_dE,97.5));

fprintf('Posterior sigma_dE: mean = %.2f, 95%% CI = [%.2f, %.2f]\n', ...
    mean(final_posterior.sig_dE), ...
    prctile(final_sig_dE,2.5), ...
    prctile(final_sig_dE,97.5));

fprintf('Posterior Ef (E0+dE): mean = %.2f, 95%% CI = [%.2f, %.2f]\n', ...
    mean(final_posterior.mean_dE + final_posterior.mean_E0),...
    prctile(final_mean_dE + final_mean_E0,2.5),...
    prctile(final_mean_dE + final_mean_E0,97.5));

fprintf('Posterior mu_tau: mean = %.2f, 95%% CI = [%.2f, %.2f]\n', ...
    mean(final_posterior.mean_tau), prctile(final_mean_tau,2.5), prctile(final_mean_tau,97.5));

fprintf('Posterior sigma_tau: mean = %.2f, 95%% CI = [%.2f, %.2f]\n', ...
    mean(final_posterior.sig_tau), prctile(final_sig_tau,2.5), prctile(final_sig_tau,97.5));

fprintf('Posterior rho: mean = %.2f, 95%% CI = [%.2f, %.2f]\n', ...
    mean(final_posterior.rho), prctile(final_rho,2.5), prctile(final_rho,97.5));

