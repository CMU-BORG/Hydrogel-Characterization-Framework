%% SMC Sampler for Cross-Sectional Area Dynamics Population
%{


%}

clear all; close all; clc

addpath(cd+"\Sequential Monte Carlo Library")
addpath(cd+"\Population Models")


% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool("local", 8);  % Adjust the number of workers if needed
end

timestamp = datestr(now, 'yyyymmdd_HHMM');

%% Loading the fitted parameter data

% data from the fitting of the exponential curves to the dynamical data.
% Used to estimate the population model parameters
data = load("AreaModelParameters.mat","params").params;
% data = data(data(:,1)==10,:);    % currently testing with just 27C data

%% Setting up the sampler and the parameters for the calibration

% Setting up the prior distributions on the statistical parameters

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

mu_A0_mean = 200; % [mm^2] (centered on the value specified in the CAD)
mu_A0_std = 50;  % [mm^2] (variance is known from printer calibration to not be too large)
mu_A0_prior = makedist("Normal","mu",mu_A0_mean, "sigma",mu_A0_std);


mu_dA_mean = 0; % [mm^2] (order of magnitude of stiffness, which for hydrogels should be ~100 kPa or 10^2 kPa)
mu_dA_std = 50;  % [kPa] 
mu_dA_prior = makedist("Normal","mu",mu_dA_mean, "sigma",mu_dA_std);

sig_A_lower = 0.00001;      % [mm^2] lowest possible value of sigma  
sig_A_upper = 100;    % [mm^2] highest possible values of sigma
sig_A_prior = makedist("Uniform","lower",sig_A_lower,"upper",sig_A_upper);

% CORRELATION COEFFICIENT PRIOR: The correlation coefficient of the
% multivariate normal is constrained to  [-1,1]. Because we don't have
% information about what it should be we can just use a uniform
% distribution prior.
rho_A_prior = makedist("Uniform", "lower", -1, "upper", 1);

% GAMMA PARAMETER PRIORS: The gamma distribution is parameterized by alpha
% and beta, but it can reparameterized by the mean and standard deviation 
% (mean = alpha*beta, variance = beta*mean). Both of these parameters have
% to be positive, so it would make sense to put gamma distribution priors
% on both. Alternatively you could use a uniform distribution on the either
% or both

% Prior for mean of tau (centered around 120 min)
mean_tau_prior = makedist("Gamma", "a", 3, "b", 40);  % mean = 120 min, std ≈ 69 min
sig_tau_prior = makedist("Gamma", "a", 2, "b", 30);  % mean = 60 min, std ≈ 42 min

% prior for the effect sizes. centered around no effect, guaranteed to be 
% within -100% change for a 20C change in temperature (both to prevent 
% negative values for variance and gamma parameters, and based on previous 
% experimentation in the lab)
effect_min = -0.9999;
effect_max = 5;
effect_prior = makedist("Uniform","lower",effect_min,"upper",effect_max); 

priors = {mu_A0_prior,mu_dA_prior,...
          sig_A_prior, sig_A_prior, rho_A_prior,...
          mean_tau_prior, sig_tau_prior,...
          effect_prior,effect_prior,effect_prior,effect_prior};

%% Step 1: Converge a starting proposal for the Normal parameters
gammas_ = linspace(0,1,10)';

% Getting an A0 distribution using all of the A0 data
model_norm_prior = GaussianModel(0);
A0_prior = SMC_Sampler(2000,model_norm_prior,priors([1,3]),2.38,data(:,2),gammas_,{"\mu_{A_0}","\sigma_{A_0}"},{},{{2,LogTransform()}},1);
A0_prior.GeneratePosterior();
A0_prior.PlotConvergence();

% Getting a dA distribution using the 27 C data
dA_27_prior = SMC_Sampler(2000,model_norm_prior,priors([2,3]),2.38,data(data(:,1)==0,3),gammas_,{"\mu_{\Delta A_{27}}","\sigma_{\Delta A_{27}}"},{},{{2,LogTransform()}},1);
dA_27_prior.GeneratePosterior();
dA_27_prior.PlotConvergence();

% Getting a dA*(1+alpha) distribution using 37 C data 
dA_37_prior = SMC_Sampler(2000,model_norm_prior,priors([2,3]),2.38,data(data(:,1)==1,3),gammas_,{"\mu_{\Delta A_{37}}","\sigma_{\Delta A_{37}}"},{},{{2,LogTransform()}},1);
dA_37_prior.GeneratePosterior();
dA_37_prior.PlotConvergence();

% estimating alpha using the estimated dA values at each temperature
alpha_mu_est = dA_37_prior.posterior(:,1) ./ dA_27_prior.posterior(:,1) - 1;
alpha_sig_est = dA_37_prior.posterior(:,2) ./ dA_27_prior.posterior(:,2) - 1;
alpha_posterior = [alpha_mu_est, alpha_sig_est];
alpha_prior = Dummy_Sampler(alpha_posterior,{"\alpha_\mu","\alpha_\sigma"},{},priors([8,9]),1);
alpha_prior.PlotPosteriorHistograms();

%% Step 2: Converge a starting proposal for the Gamma parameters

% Getting a tau distribution using the 27 C data
model_gamma_prior = GammaModel(1);

% estimating tau parameters for 27 C data
tau_27_prior = SMC_Sampler(2000,model_gamma_prior,priors([6,7]),2.38,data(data(:,1)==0,4),gammas_,{"\alpha_{\mu_{27}}","\alpha_{\sigma_{27}}"},{},{{1,LogTransform()},{2,LogTransform()}},1);
tau_27_prior.GeneratePosterior();
tau_27_prior.PlotConvergence();

% estimating tau parameters for 37 C data
tau_37_prior = SMC_Sampler(2000,model_gamma_prior,priors([6,7]),2.38,data(data(:,1)==1,4),gammas_,{"\alpha_{\mu_{37}}","\alpha_{\sigma_{37}}"},{},{{1,LogTransform()},{2,LogTransform()}},1);
tau_37_prior.GeneratePosterior();
tau_37_prior.PlotConvergence();

% estimating beta using the estimate tau at each temperature
beta_mu_est = tau_37_prior.posterior(:,1) ./ tau_27_prior.posterior(:,1) - 1;
beta_sig_est = tau_37_prior.posterior(:,2) ./ tau_27_prior.posterior(:,2) - 1;
beta_posterior = [beta_mu_est, beta_sig_est];
beta_prior = Dummy_Sampler(beta_posterior,{"\beta_\mu","\beta_\sigma"},{},priors([10,11]),1);

beta_prior.PlotPosteriorHistograms();

%% Step 3: Construct the sampler using the pre-convergence proposal distributions to get full joint posterior
gammas_ = linspace(0.75,1,5)'; % because we already have a good starting place, we don't have to start all the way at the prior again, but can start most of the way to the posterior

parameter_labels = {"\mu_{A_0}","\mu_{\Delta A}",...
                    "\sigma_{A_0}","\sigma_{\Delta A}","\rho",...
                    "\mu_\tau","\sigma_\tau",...
                    "\alpha_\mu","\alpha_\sigma",...
                    "\beta_\mu","\beta_\sigma"};

model = PopulationModel(1); % ~~~~ NEED TO CHANGE TO 1 ONCE THE 37C DATA IS INCLUDED ~~~~~
mappings = {1,[1,3]; % A0 sampler 
            2,[2,4]; % dA_27 sampler
            3,[5];   % rho prior
            4,[6,7]; % tau 27 sampler
            5,[8,9]; % alpha sampler
            6,[10,11]}; % beta sampler
priors = {A0_prior, dA_27_prior, rho_A_prior, tau_27_prior, alpha_prior, beta_prior};

transforms = {{3, LogTransform()};
              {4, LogTransform()};
              {5, LogitTransform(-1,1)};
              {6, LogTransform()};
              {7, LogTransform()};
              {8, LogitTransform(effect_min,effect_max)};
              {9, LogitTransform(effect_min,effect_max)};
              {10, LogitTransform(effect_min,effect_max)};
              {11, LogitTransform(effect_min,effect_max)};};

sampler = SMC_Sampler(5000,model,priors,2.38,data,gammas_,parameter_labels,mappings,transforms,0);
sampler.GeneratePosterior();

sampler.PlotConvergence();

%%
summary = sampler.PosteriorSummary();

%%
sampler.PlotPosteriorHistograms();
% sampler.PlotJointPosterior();

%% Mean Posterior Prediction
theta_post = summary(:,1);

model.setParams(theta_post);

t = 0:15:360;
N_traj = 100000;

figure("Position",[100,100,600,400],"Color","w"); hold all

[trajectories,particles] = model.rnd_trajectories(t,N_traj,0);
mean_27 = mean(trajectories,1);
std_27 = std(trajectories,[],1);
% for i=1:N_traj
%     plot(t,trajectories(i,:),"Color",[1,0.5,0.5]*rand(1),"LineWidth",0.1)
% end
errorbar(t,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','r','Color','r')

[trajectories,particles] = model.rnd_trajectories(t,N_traj,1);
mean_37 = mean(trajectories,1);
std_37 = std(trajectories,[],1);
% for i=1:N_traj
%     plot(t,trajectories(i,:),"Color",[1,0.5,0.5]*rand(1),"LineWidth",0.1)
% end
errorbar(t,mean_37,std_37,std_37,'vertical','Marker','s','MarkerFaceColor','b','Color','b')
