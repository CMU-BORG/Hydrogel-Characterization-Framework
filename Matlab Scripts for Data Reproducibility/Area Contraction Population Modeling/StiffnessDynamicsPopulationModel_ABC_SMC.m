%% SMC Sampler for Stiffness Dynamics Population
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

%% Loading the Experimental data

% data location
% currentFolder = pwd;
root = "../Experimental Data/";

% loading the inventory table with all of the summary data
inventory_file = "Inventory Final.xlsx";
inventory = readtable(root + inventory_file,"ReadVariableNames",true);

% pulling out relevant columns of the inventory table (CrosslinkTime,
% CrosslinkTemp, EffectiveElasticModulus)
data = inventory(:, [4, 3, 16]);

rows27 = find(inventory.CrosslinkTemp == 27)';
rows37 = find(inventory.CrosslinkTemp == 37)';
rows0 = find(inventory.CrosslinkTemp == 0)'; % for samples that were tested immediately with no additional cross-linking

rows = [rows0, rows27, rows37];
data = data{rows,:}; % 27C and 37C Initial and Final mass data
data(:,2) = data(:,2)*60;   % converting to minutes for agreement with the area data units
data(:,3) = data(:,3)*1000; % converting to kPa

data(data(:,1)==0,1) = -1; 
data(data(:,1)>0,1) = data(data(:,1)>0,1)==37;
data_27 = data((data(:,1)==0)|(data(:,1)==-1),:);
data_37 = data((data(:,1)==1)|(data(:,1)==-1),:);
%%
figure; hold all
% plot(0*data(:,2),data(:,3),'^k')
% plot(data_27(:,2),data_27(:,4),'sr');
% plot(data_37(:,2),data_37(:,4),'ob');
% errorbar(0*data(:,2),median(data(:,3)), quantile(data(:,3),[0.25]) - median(data(:,3)), quantile(data(:,3),[0.75]) - median(data(:,3)),'sr',"LineWidth",2)

for i=60*[0,2,4,6]
    inds = find(data_27(:,2)==i);
    data_i = data_27(inds,3);
    errorbar(data_27(inds,2),mean(data_i), std(data_i), std(data_i),'sr',"LineWidth",4)
    % errorbar(data_27(inds,2),median(data_i), quantile(data_i,[0.025]) - median(data_i), quantile(data_i,[0.975]) - median(data_i),'k',"LineWidth",2,"CapSize",0)
    
    inds = find(data_37(:,2)==i);
    data_i = data_37(inds,3);
    errorbar(data_37(inds,2),mean(data_i), std(data_i), std(data_i),'sb',"LineWidth",4)
    % errorbar(data_37(inds,2),median(data_i), quantile(data_i,[0.25]) - median(data_i), quantile(data_i,[0.75]) - median(data_i),'sb',"LineWidth",2)
    
end
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

% prior for the effect sizes. centered around no effect, guaranteed to be 
% within -100% change for a 20C change in temperature (both to prevent 
% negative values for variance and gamma parameters, and based on previous 
% experimentation in the lab)
effect_min = -0.9999;
effect_max = 5;
effect_prior = makedist("Uniform","lower",effect_min,"upper",effect_max); 

priors = {mu_E0_prior,mu_dE_prior,...
          sig_E_prior, sig_E_prior, rho_E_prior,...
          mean_tau_prior, sig_tau_prior,...
          effect_prior,effect_prior,effect_prior,effect_prior};

parameter_labels = {"\mu_{E_0}","\mu_{\Delta E}",...
                    "\sigma_{E_0}","\sigma_{\Delta E}","\rho",...
                    "\mu_\tau","\sigma_\tau",...
                    "\alpha_\mu","\alpha_\sigma",...
                    "\beta_\mu","\beta_\sigma"};


<<<<<<< Updated upstream
%% Step 1: Estimate the initial distribution using standard SMC
gammas_ = (linspace(0,1,5)');

transforms = cell(1,1);
transforms{1} = {2, LogTransform()};

% Getting an E0 distribution using all of the E0 data
model_norm_prior = GaussianModel(0);
E0_prior = SMC_Sampler(2000,model_norm_prior,priors([1,3]),2.38,data(data(:,1)==-1,3),gammas_,{"\mu_{E_0}","\sigma_{E_0}"},{},transforms,1);
E0_prior.GeneratePosterior();
E0_prior.PlotConvergence();

summary = E0_prior.PosteriorSummary();
=======
%% Step 1: Estimate the initial mass distribution using standard SMC
% gammas_ = (linspace(0,1,5)');
% 
% transforms = cell(1,1);
% transforms{1} = {2, LogTransform()};
% 
% % Getting an M0 distribution using all of the M0 data
% model_norm_prior = GaussianModel(0);
% E0_prior = SMC_Sampler(2000,model_norm_prior,priors([1,3]),2.38,data(data(:,1)==-1,3),gammas_,{"\mu_{E_0}","\sigma_{E_0}"},{},transforms,1);
% E0_prior.GeneratePosterior();
% E0_prior.PlotConvergence();
% 
% summary = E0_prior.PosteriorSummary();

>>>>>>> Stashed changes

%% Step 2: Estimating the Population model for 27C and 37C separately
%{
    Don't spend too much time here. Only get down to 7.5% error on the
    means and 15% error on the standard deviations. Then with the
    partially-calibrated sampling distribution, do the full run on the
    whole model and data set, pushing the error bounds more.
%}

stiff_pop_model_one_temp = Stiffness_PopulationModel(0);

N_gens = 10;
gammas_ = linspace(0,1,N_gens)';

% 27 C data

% constructing priors and mappings
priors_ = priors(1:7); 
mappings_ = {};
transforms = {{3,LogTransform()}; {4,LogTransform()}; {5,LogitTransform(-1,1)}; {6, LogTransform()}; {7, LogTransform()}};

data_27_prior = SMC_Sampler(2000, stiff_pop_model_one_temp, priors_,2,data_27,gammas_,parameter_labels(1:7), mappings_, transforms,1);

data_27_prior.GeneratePosterior();
data_27_prior.PlotConvergence();

fprintf("\t~~~~ 27 C Data ~~~~\n")
summary = data_27_prior.PosteriorSummary();

%% Posterior Predictive Check 27C Data
theta_post = summary(:,1);

stiff_pop_model_one_temp.setParams(theta_post);

t = 0:15:360;
N_traj = 100000;

figure("Position",[100,100,600,400],"Color","w"); hold all

[trajectories,particles] = stiff_pop_model_one_temp.rnd_trajectories(t,N_traj,0);

plot(data_27(:,2),data_27(:,3),'ok')
t_groups = unique(data_27(:,2));
N_samps_per_group = 6;
for i=1:length(t_groups)
    inds = data_27(:,2)==t_groups(i);
    data_ = data_27(inds,3);
    errorbar(t_groups(i),mean(data_), std(data_), std(data_),'sg',"LineWidth",2)
end

mean_27 = mean(trajectories,1);
std_27 = std(trajectories,[],1);
errorbar(t,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','r','Color','r')


%% 37 C data

data_37_prior = SMC_Sampler(2000, stiff_pop_model_one_temp, priors_,2,data_37,gammas_,parameter_labels(1:7), mappings_, transforms,1);
data_37_prior.GeneratePosterior();
data_37_prior.PlotConvergence();

fprintf("\t~~~~ 37 C Data ~~~~\n")
summary = data_37_prior.PosteriorSummary();

%% Posterior Predictive Check 37C Data
theta_post = summary(:,1);

stiff_pop_model_one_temp.setParams(theta_post);

t = 0:15:360;
N_traj = 100000;

figure("Position",[100,100,600,400],"Color","w"); hold all

[trajectories,particles] = stiff_pop_model_one_temp.rnd_trajectories(t,N_traj,0);

plot(data_37(:,2),data_37(:,3),'ok')
t_groups = unique(data_37(:,2));
N_samps_per_group = 6;
mean_data = zeros(length(t_groups),1);
for i=1:length(t_groups)
    inds = data_37(:,2)==t_groups(i);
    data_ = data_37(inds,3);
    mean_data(i) = mean(data_);
    errorbar(t_groups(i),mean(data_), std(data_), std(data_),'sg',"LineWidth",2)
end

mean_37 = mean(trajectories,1);
std_37 = std(trajectories,[],1);
errorbar(t,mean_37,std_37,std_37,'vertical','Marker','s','MarkerFaceColor','b','Color','b')

%%
% estimating alpha and beta using the estimated dA values at each temperature
alpha_mu_est = min(max(data_37_prior.posterior(:,2) ./ data_27_prior.posterior(:,2) - 1, effect_min), 0.99999*effect_max);
alpha_sig_est = min(max(data_37_prior.posterior(:,4) ./ data_27_prior.posterior(:,4) - 1, effect_min), 0.99999*effect_max);
beta_mu_est = min(max(data_37_prior.posterior(:,6) ./ data_27_prior.posterior(:,6) - 1, effect_min), 0.99999*effect_max);
beta_sig_est = min(max(data_37_prior.posterior(:,7) ./ data_27_prior.posterior(:,7) - 1, effect_min), 0.99999*effect_max);

effect_posterior = [alpha_mu_est, alpha_sig_est, beta_mu_est, beta_sig_est];
effect_prior = Dummy_Sampler(effect_posterior,{"\alpha_\mu","\alpha_\sigma","\beta_\mu","\beta_\sigma"},{},priors([8,9,10,11]),1);
effect_prior.PlotPosteriorHistograms();

%% Step 3: Estimate the full model

stiffness_model = Stiffness_PopulationModel(1);

mappings = {1,[1,2,3,4,5,6,7];    % initial parameter estimates
            2,[8,9,10,11]};     % effect parameter estimates

priors_final = {data_27_prior, effect_prior};

transforms = {{3, LogTransform()};
              {4, LogTransform()};
              {5, LogitTransform(-1,1)};
              {6, LogTransform()};
              {7, LogTransform()};
              {8, LogitTransform(-1,effect_max)};
              {9, LogitTransform(-1,effect_max)};
              {10, LogitTransform(-1,effect_max)};
              {11, LogitTransform(-1,effect_max)};};

% all data
% data(data(:,1)==-1,1) = 0;
% 
% temp_vals = unique(data(:,1));
% data_cell = cell(7,1);
% for j=1:length(temp_vals)
%     data_ = data(data(:,1)==temp_vals(j),:);
%     time_vals = unique(data_(:,2));
% 
%     for i=1:length(time_vals)
%         data_i = data_(data_(:,2)==time_vals(i),:);
%         data_cell{(j-1)*4 + i} = {length(data_i), [time_vals(i),temp_vals(j)], data_i(:,3)};
%     end
% end
% 
% N_sims = 5;
% gammas_mu = linspace(0.3,0.1,N_gens)';
% gammas = zeros(N_gens,1,5);
% var_ind = [5];
% for i=1:5
%     gammas(:,1,i) = gammas_mu*(0.5 + 0.5*(any(i==var_ind)));
% end
gammas_ = linspace(0.5,1,10)';

%%
% sampler = ABC_SMC_Sampler(4000,N_sims,stiffness_model,priors_final,2.38,data,gammas,parameter_labels,mappings,transforms,0);
sampler = SMC_Sampler(2000, stiffness_model, priors_final,1.5,data,gammas_,parameter_labels, mappings, transforms,0);

sampler.GeneratePosterior();

summary = sampler.PosteriorSummary();

sampler.PlotConvergence();

sampler.PlotPosteriorHistograms();

%% Mean Posterior Prediction
theta_post = summary(:,1);

stiffness_model.setParams(theta_post);

t = 0:15:360;
N_traj = 100000;

figure("Position",[100,100,600,400],"Color","w"); hold all

[trajectories,particles] = stiffness_model.rnd_trajectories(t,N_traj,0);
mean_27 = mean(trajectories,1);
std_27 = std(trajectories,[],1);
% for i=1:N_traj
%     plot(t,trajectories(i,:),"Color",[1,0.5,0.5]*rand(1),"LineWidth",0.1)
% end
errorbar(t,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','r','Color','r')

[trajectories,particles] = stiffness_model.rnd_trajectories(t,N_traj,1);
mean_37 = mean(trajectories,1);
std_37 = std(trajectories,[],1);
% for i=1:N_traj
%     plot(t,trajectories(i,:),"Color",[1,0.5,0.5]*rand(1),"LineWidth",0.1)
% end
errorbar(t,mean_37,std_37,std_37,'vertical','Marker','s','MarkerFaceColor','b','Color','b')
