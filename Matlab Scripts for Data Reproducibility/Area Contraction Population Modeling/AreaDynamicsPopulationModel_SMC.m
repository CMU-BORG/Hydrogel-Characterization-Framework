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
% sampler.PlotPosteriorHistograms();
sampler.PlotKDEstimate()
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


%% Compare first and last generation particle trajectories for one temperature

tempCondition = 0;  % 0 = 27C, 1 = 37C
dataTemp = data(data(:,1) == tempCondition, :);
tPlot = linspace(0, 360, 200);
nShow = 100;
rng(1);  % reproducible selection

% particleHistory = sampler.particle_history;
% 
% firstGenParticles = squeeze(particleHistory(1, :, :));   % 5000 x 11
% lastGenParticles  = squeeze(particleHistory(end, :, :)); % 5000 x 11

% Build a synthetic pre-convergence ensemble from stage-1/2 posterior samples
nPrePool = size(A0_prior.posterior, 1);
preGenParticles = zeros(nPrePool, 11);

for k = 1:nPrePool
    iA = randi(size(A0_prior.posterior, 1));
    iD = randi(size(dA_27_prior.posterior, 1));
    iT = randi(size(tau_27_prior.posterior, 1));
    iAlpha = randi(size(alpha_prior.posterior, 1));
    iBeta = randi(size(beta_prior.posterior, 1));

    preGenParticles(k,:) = [
        A0_prior.posterior(iA,1), ...
        dA_27_prior.posterior(iD,1), ...
        A0_prior.posterior(iA,2), ...
        dA_27_prior.posterior(iD,2), ...
        random(rho_A_prior), ...
        tau_27_prior.posterior(iT,1), ...
        tau_27_prior.posterior(iT,2), ...
        alpha_prior.posterior(iAlpha,1), ...
        alpha_prior.posterior(iAlpha,2), ...
        beta_prior.posterior(iBeta,1), ...
        beta_prior.posterior(iBeta,2)];
end

% Final post-convergence generation from the joint sampler
postGenParticles = squeeze(sampler.particle_history(end, :, :));

idxPre = randperm(size(preGenParticles, 1), min(nShow, size(preGenParticles, 1)));
idxPost = randperm(size(postGenParticles, 1), min(nShow, size(postGenParticles, 1)));
% idxFirst = randperm(size(firstGenParticles, 1), min(nShow, size(firstGenParticles, 1)));
% idxLast  = randperm(size(lastGenParticles, 1),  min(nShow, size(lastGenParticles, 1)));

figure("Color","w","Position",[100,100,700,450]); hold on

% plot first gen (pre-conv)

% nFirst = numel(idxFirst);
% endValsFirst = zeros(nFirst, 1);
% trajStoreFirst = zeros(nFirst, numel(tPlot));
nPre = numel(idxPre);
endValsPre = zeros(nPre, 1);
trajStorePre = zeros(nPre, numel(tPlot));

for k = 1:nPre %nFirst
    % theta = firstGenParticles(idxFirst(k), :).';
    % model.setParams(theta);
    % [traj, ~] = model.rnd_trajectories(tPlot, 1, tempCondition);
    % 
    % trajNorm = traj(1,:) ./ traj(1,1);
    % endValsFirst(k) = trajNorm(end);
    % trajStoreFirst(k,:) = trajNorm;

    theta = preGenParticles(idxPre(k), :).';
    model.setParams(theta);
    [traj, ~] = model.rnd_trajectories(tPlot, 1, tempCondition);

    trajNorm = traj(1,:) ./ traj(1,1);
    endValsPre(k) = trajNorm(end);
    trajStorePre(k,:) = trajNorm;
end

% qFirst = quantile(endValsFirst, [0.25 0.75]);
% keepFirst = endValsFirst >= qFirst(1) & endValsFirst <= qFirst(2);
% plot(tPlot, trajStoreFirst(keepFirst,:)', "Color", [0.82 0.82 0.82], "LineWidth", 0.8);

qPre = quantile(endValsPre, [0.25 0.75]);
keepPre = endValsPre >= qPre(1) & endValsPre <= qPre(2);
plot(tPlot, trajStorePre(keepPre,:)', "Color", [0.82 0.82 0.82], "LineWidth", 0.8);

% plot last gen (post-conv)
nPost = numel(idxPost);
endValsPost = zeros(nPost, 1);
trajStorePost = zeros(nPost, numel(tPlot));
% nLast = numel(idxLast);
% endValsLast = zeros(nLast, 1);
% trajStoreLast = zeros(nLast, numel(tPlot));

for k = 1:nPost
    % theta = lastGenParticles(idxLast(k), :).';
    % model.setParams(theta);
    % [traj, ~] = model.rnd_trajectories(tPlot, 1, tempCondition);
    % 
    % trajNorm = traj(1,:) ./ traj(1,1);
    % endValsLast(k) = trajNorm(end);
    % trajStoreLast(k,:) = trajNorm;

    theta = postGenParticles(idxPost(k), :).';
    model.setParams(theta);
    [traj, ~] = model.rnd_trajectories(tPlot, 1, tempCondition);

    trajNorm = traj(1,:) ./ traj(1,1);
    endValsPost(k) = trajNorm(end);
    trajStorePost(k,:) = trajNorm;
end

% qLast = quantile(endValsLast, [0.25 0.75]);
% keepLast = endValsLast >= qLast(1) & endValsLast <= qLast(2);
% plot(tPlot, trajStoreLast(keepLast,:)', "Color", [0.10 0.35 0.80], "LineWidth", 0.8);
% 
% model.setParams(theta_post);
% [fitTraj, ~] = model.rnd_trajectories(tPlot, 1, tempCondition);
% fitNorm = fitTraj(1,:) ./ fitTraj(1,1);
% plot(tPlot, fitNorm, "k-", "LineWidth", 2.5);

qPost = quantile(endValsPost, [0.25 0.75]);
keepPost = endValsPost >= qPost(1) & endValsPost <= qPost(2);
plot(tPlot, trajStorePost(keepPost,:)', "Color", [0.10 0.35 0.80], "LineWidth", 0.8);

% % Mean posterior fit from the final joint sampler
% summary = sampler.PosteriorSummary();
% theta_post = summary(:,1);
% 
% model.setParams(theta_post);
% [fitTraj, ~] = model.rnd_trajectories(tPlot, 1, tempCondition);
% fitNorm = fitTraj(1,:) ./ fitTraj(1,1);
% plot(tPlot, fitNorm, "k-", "LineWidth", 2.5);

fitNorm = median(trajStorePost, 1);
plot(tPlot, fitNorm, "k-", "LineWidth", 2.5);
xlim([0, 360]);

xlabel("Crosslink time [min]")
ylabel("Normalized Cross-sectional Area")
title("First and last generation particle trajectories for 27C")
box on
legend({"IQR of First Gen", "IQR of Last Gen", "Median Post-Gen Trajectory"}, "Location", "best")


%% Save all figures
% resultsDir = fullfile(pwd, "SMC_results_" + string(timestamp));
% if ~exist(resultsDir, "dir")
%     mkdir(resultsDir);
% end
% 
% figNames = {
%     "A0_mu_convergence.svg"
%     "A0_sigma_convergence.svg"
%     "dA27_mu_convergence.svg"
%     "dA27_sigma_convergence.svg"
%     "dA37_mu_convergence.svg"
%     "dA37_sigma_convergence.svg"
%     "alpha_mu_histogram.svg"
%     "alpha_sigma_histogram.svg"
%     "tau27_mu_convergence.svg"
%     "tau27_sigma_convergence.svg"
%     "tau37_mu_convergence.svg"
%     "tau37_sigma_convergence.svg"
%     "beta_mu_histogram.svg"
%     "beta_sigma_histogram.svg"
%     "full_mu_A0_convergence.svg"
%     "full_mu_dA_convergence.svg"
%     "full_sigma_A0_convergence.svg"
%     "full_sigma_dA_convergence.svg"
%     "full_rho_convergence.svg"
%     "full_mu_tau_convergence.svg"
%     "full_sigma_tau_convergence.svg"
%     "full_alpha_mu_convergence.svg"
%     "full_alpha_sigma_convergence.svg"
%     "full_beta_mu_convergence.svg"
%     "full_beta_sigma_convergence.svg"
%     "full_mu_A0_histogram.svg"
%     "full_mu_dA_histogram.svg"
%     "full_sigma_A0_histogram.svg"
%     "full_sigma_dA_histogram.svg"
%     "full_rho_histogram.svg"
%     "full_mu_tau_histogram.svg"
%     "full_sigma_tau_histogram.svg"
%     "full_alpha_mu_histogram.svg"
%     "full_alpha_sigma_histogram.svg"
%     "full_beta_mu_histogram.svg"
%     "full_beta_sigma_histogram.svg"
%     "mean_posterior_prediction.svg"
%     "particle_generation_comparison_27C.svg"
%     };
% 
% figs = findall(groot, "Type", "figure");
% [~, idx] = sort(arrayfun(@(f) f.Number, figs));
% figs = figs(idx);
% 
% if numel(figs) ~= numel(figNames)
%     error("Found %d figures, but expected %d.", numel(figs), numel(figNames));
% end
% 
% for k = 1:numel(figs)
%     exportgraphics(figs(k), fullfile(resultsDir, figNames{k}), "ContentType", "vector");
% end
