%% Examining Summary Statistics for the Population Model as a function of parameter
%{

%}

close all; clear all; clc

model = Stiffness_PopulationModel(0);

%% Setting up the priors

mu_E0_mean = 100; % [kPa] (around typical values for hydrogel)
mu_E0_std = 50;  % [kpa] (broad prior)
mu_E0_prior = makedist("Normal","mu",mu_E0_mean, "sigma",mu_E0_std);

% we know the stiffness increases with cross-linking but we are uninformed
% about the magnitude
mu_dE_std = 200;  % [g] 
mu_dE_prior = makedist("HalfNormal", "sigma",mu_dE_std);

% sig_E_lower = 1;      % [g] lowest possible value of sigma  
% sig_E_upper = 100;    % [g] highest possible values of sigma
% sig_E_prior = makedist("Uniform","lower",sig_E_lower,"upper",sig_E_upper);
sig_E_prior = makedist("HalfNormal",'sigma',50);

% CORRELATION COEFFICIENT PRIOR: 
rho_E_prior = makedist("Uniform", "lower", -1, "upper", 1);

% GAMMA PARAMETER PRIORS: The gamma distribution is parameterized by alpha
% and beta, but it can reparameterized by the mean and standard deviation 
% (mean = alpha*beta, variance = beta*mean). Both of these parameters have
% to be positive, so it would make sense to put gamma distribution priors
% on both. Alternatively you could use a uniform distribution on the either
% or both

% Prior for mean of tau (centered around 120 min)
mean_tau_prior = makedist("Gamma", "a", 3, "b", 40);  % mean = 120 min, std ≈ 69 min
sig_tau_prior = makedist("Gamma", "a", 2, "b", 30);  % mean = 60 min, std ≈ 42 min


priors = {mu_E0_prior,mu_dE_prior,...
          sig_E_prior, sig_E_prior, rho_E_prior,...
          mean_tau_prior, sig_tau_prior};

parameter_labels = {"\mu_{E_0}","\mu_{\Delta E}",...
                    "\sigma_{E_0}","\sigma_{\Delta E}","\rho",...
                    "\mu_\tau","\sigma_\tau"};

prior = Dummy_Sampler(zeros(1000,7),parameter_labels,{},priors,1);

N_parameter_samples = 10000;
N_data_resamples = 10;

% test run to get the size of the summary stats
theta_i = prior.SamplePrior(1);
data_j = RunExperiment(model,theta_i);
sm_sts_j = CalculateSummaryStatistics(data_j);


thetas_ = zeros(N_data_resamples * N_parameter_samples, 7);
sm_sts_ = zeros(N_data_resamples * N_data_resamples, length(sm_sts_j));

index = 1;
for i=1:N_parameter_samples
    theta_i = prior.SamplePrior(1);
    for j=1:N_data_resamples
        data_j = RunExperiment(model,theta_i);
        sm_sts_j = CalculateSummaryStatistics(data_j);
        thetas_(index,:) = theta_i;
        sm_sts_(index,:) = sm_sts_j;
        index = index + 1;
    end
end
   
% Multi-Linear Regression on the simulated summary statistics
ONE = ones(N_data_resamples * N_parameter_samples,1);
A = [sm_sts_];
AA = A'*A; I = eye(size(AA));
beta_ = (AA + 0.1*I)\(A'*thetas_);

theta_est = @(sm_sts) [sm_sts] * beta_ ;

%% Summarizing the regression
theta_est_all = theta_est(sm_sts_);
error_all = 100*abs(theta_est_all - thetas_)./abs(thetas_);
median_error = median(error_all,1)

%%
% theta_i = prior.SamplePrior(1)
% data_j = RunExperiment(model,theta_i);
% sm_sts_j = CalculateSummaryStatistics(data_j);
% est_theta_i = theta_est(sm_sts_j')
% err = 100*(theta_i - est_theta_i)./theta_i


function data = RunExperiment(model,theta)
%{
    Using the provided model and parameter set theta, run a simulated
    experiment to provide the same data format as in the main study.
%}
    model.setParams(theta);

    t = 0:120:360;
    N_traj = 6;
    data = zeros(N_traj,length(t));
    
    for i=1:length(t)
        [data(:,i),~] = model.rnd_trajectories(t(i),N_traj,0);
    end
end

function summary_stats = CalculateSummaryStatistics(data)
%{
    Based on set of data from a simulated experiment, calculate a large
    vector of possible summary statistics that are then regressed to
    estimated the true values used to generate the data in the first place.
%}

summary_stats = zeros(18,1);

index = 1;
N_time_points = 4;

%% time point level summaries
% time-point-wise means
for i=1:N_time_points
    summary_stats(index) = mean(data(:,i));
    index = index+1;
end

% time-point-wise standard deviations
for i=1:N_time_points
    summary_stats(index) = std(data(:,i));
    index = index+1;
end

% time-point-wise minimum
for i=1:N_time_points
    summary_stats(index) = min(data(:,i));
    index = index+1;
end

% time-point-wise maximum
for i=1:N_time_points
    summary_stats(index) = max(data(:,i));
    index = index+1;
end

%% Global level summaries
summary_stats(index) = median(data,"all");
index = index + 1;
summary_stats(index) = std(data,[],"all");
index = index + 1;

%% What is we just used the data itself
summary_stats = [summary_stats; data(:); data(:).^2; data(:).^3; data(:).^4];

%% Adding in nonlinear elements
% summary_stats = [summary_stats; summary_stats.^2; summary_stats.^3; summary_stats.^4];

end