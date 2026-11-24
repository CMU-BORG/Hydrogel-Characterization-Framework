%% Posterior Inference Using Calibrated Model
%{
    Perform median posterior model inference regarding time to steady
    state.

    Last Updated: 11/24/2026 MJB
%}


clear all; close all; clc

addpath(cd+"\Sequential Monte Carlo Library")
addpath(cd+"\Population Models")

sampler_file_name = "SMC_sampler_results.mat";

%% CSA Data
area_folder = "Area Results20261108_0000\";
area_data = load(area_folder + sampler_file_name).sampler;

% using the calibrated median parameters construct the population model for
% the CSA dynamics
posterior_median = area_data.PosteriorSummary();
area_model = PopulationModel(1);
area_model.setParams(posterior_median(:,1));

% generate samples for the low temp group (27 C)
samps_low_temp = area_model.rnd_params(100000,0);
A0_samps = samps_low_temp(:,1);
taus_low_temp = samps_low_temp(:,3);

% generate samples for the high temp group (37 C)
samps_high_temp = area_model.rnd_params(100000,1);
taus_high_temp = samps_high_temp(:,3);

% Calculate the time to reach steady state for each group
steady_state_frac = 0.95;
steady_state_time_low = -taus_low_temp.*log(1-steady_state_frac) / 60;
steady_state_time_high = -taus_high_temp.*log(1-steady_state_frac) / 60;

figure; hold all
histogram(steady_state_time_low,'Normalization','pdf','NumBins',50)
histogram(steady_state_time_high,'Normalization','pdf','NumBins',50)
xlabel("Time to Reach Steady State [hrs]")
legend("27 C","37 C")

fprintf("Time to Reach Steady State [hrs]:\n")
fprintf("\t27C: %.2f [%.2f, %.2f]\n",mean(steady_state_time_low),quantile(steady_state_time_low,0.025),quantile(steady_state_time_low,0.975))
fprintf("\t37C: %.2f [%.2f, %.2f]\n",mean(steady_state_time_high),quantile(steady_state_time_high,0.025),quantile(steady_state_time_high,0.975))