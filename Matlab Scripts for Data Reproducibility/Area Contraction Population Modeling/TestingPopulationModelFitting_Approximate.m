%% Testing the Synthetic Likelihood Approach to Fitting the Population Model
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

%% Generating Synthetic Data using the Population Model

N_samps_per_group = 6; % use a lot of data to make sure its right
t_groups = 0:120:360;

data = zeros(N_samps_per_group*length(t_groups), 3);

model_true = Stiffness_PopulationModel(0);

% rng(1000000); % setting for reproducibility across runs
mu_E0_true = 100; sig_E0_true = 20;
mu_dE_true = 50;  sig_dE_true = 15;
rho_true = -0.;
mu_tau_true = 150; sig_tau_true = 20;
theta_true = [mu_E0_true,mu_dE_true,sig_E0_true,sig_dE_true,rho_true,mu_tau_true,sig_tau_true];

model_true.setParams(theta_true);

for i=1:length(t_groups)
    inds = 1+N_samps_per_group*(i-1):N_samps_per_group*i;
    data(inds,1) = -1*(t_groups(i)==0);
    data(inds,2) = t_groups(i);
    data(inds,3) = model_true.rnd_trajectories(t_groups(i),N_samps_per_group,0);
end

figure(); hold all
plot(data(:,2),data(:,3),'ok')
for i=1:length(t_groups)
    inds = 1+N_samps_per_group*(i-1):N_samps_per_group*i;
    data_ = data(inds,3);
    errorbar(t_groups(i),median(data_), -std(data_), std(data_),'sr',"LineWidth",2)
end

%% Setting up the priors
mu_E0_mean = 100; % [kPa] (around typical values for hydrogel)
mu_E0_std = 50;  % [kpa] (broad prior)
mu_E0_prior = makedist("Normal","mu",mu_E0_mean, "sigma",mu_E0_std);

% we know the stiffness increases with cross-linking but we are uninformed
% about the magnitude
mu_dE_std = 100;  % [g] 
mu_dE_prior = makedist("HalfNormal", "sigma",mu_dE_std);

sig_E_prior = makedist("HalfNormal",'sigma',50);

% CORRELATION COEFFICIENT PRIOR: 
rho_E_prior = makedist("Uniform", "lower", -1, "upper", 1);

% Prior for mean of tau (centered around 120 min)
mean_tau_prior = makedist("Gamma", "a", 3, "b", 40);  % mean = 120 min, std ≈ 69 min
sig_tau_prior = makedist("Gamma", "a", 2, "b", 30);  % mean = 60 min, std ≈ 42 min

priors = {mu_E0_prior,mu_dE_prior,...
          sig_E_prior, sig_E_prior, rho_E_prior,...
          mean_tau_prior, sig_tau_prior};

parameter_labels = {"\mu_{E_0}","\mu_{\Delta E}",...
                    "\sigma_{E_0}","\sigma_{\Delta E}","\rho",...
                    "\mu_\tau","\sigma_\tau"};

%% Fitting the Model Parameters

model = Stiffness_PopulationModel(0);

N_gens = 10;
gammas_ = linspace(0,1,N_gens)';

mappings_ = {};
transforms = {{2,LogTransform()}; {3,LogTransform()}; {4,LogTransform()}; {5,LogitTransform(-1,1)}; {6, LogTransform()}; {7, LogTransform()}};

sampler = SMC_Sampler(2000, model, priors,2,data,gammas_,parameter_labels, mappings_, transforms,0);

sampler.GeneratePosterior();

%% Histograms
sampler.PlotPosteriorHistograms()

%% Plotting the Results

figs = sampler.PlotConvergence();
for i=1:length(figs)
    figure(figs{i}); hold all
    yline(theta_true(i),'--b','LineWidth',3)
end
fprintf("\t~~~~ Estimated Parameters ~~~~\n")
summary = sampler.PosteriorSummary();

%% Posterior Predictive Check
theta_post = summary(:,1);

t = 0:15:360;
N_traj = 100000;

model.setParams(theta_post);
[trajectories,particles] = model.rnd_trajectories(t,N_traj,0);
model.setParams(theta_true);
[trajectories_true,particles] = model.rnd_trajectories(t,N_traj,0);

figure("Position",[100,100,600,400],"Color","w"); hold all

plot(data(:,2),data(:,3),'ok')
for i=1:length(t_groups)
    inds = 1+N_samps_per_group*(i-1):N_samps_per_group*i;
    data_ = data(inds,3);
    errorbar(t_groups(i),median(data_), -std(data_), std(data_),'sg',"LineWidth",2)
end

mean_27 = mean(trajectories,1);
std_27 = std(trajectories,[],1);
errorbar(t,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','r','Color','r')

mean_27 = mean(trajectories_true,1);
std_27 = std(trajectories_true,[],1);
errorbar(t,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','b','Color','b')


%% Testing the Assumption of Mean Dynamics

% ft = @(x,t_) x(1).*(1 - exp(-t_./x(2))); 
% t_ = linspace(0,400,100);
% 
% replicates = 1;
% bias = zeros(replicates,3);
% for j=1:replicates
%     t = 0:120:360;
%     N_traj = 100;
% 
%     particles = zeros(N_traj*length(t),3);
%     trajectories_true = zeros(N_traj,length(t));
% 
%     for i=1:length(t)
%         [trajectories_true(:,i),particles(1+(i-1)*N_traj:i*N_traj,:)] = model.rnd_trajectories(t(i),N_traj,0);
%     end
% 
%     mean_trajectory = mean(trajectories_true,1);
%     E0 = mean_trajectory(1);
%     [dE, tau] = model.FitRelaxation(t,mean_trajectory' - E0);
%     fit_params = [E0,dE,tau];
%     true_params = [mean(particles(:,1)),mean(particles(:,2)),1/mean(1./particles(:,3))];
%     bias_j = 100*(true_params - fit_params)./true_params;
%     bias(j,:) = bias_j;
%     figure; hold all
%     plot(t,trajectories_true,'ok')
%     plot(t,mean_trajectory,'-r','LineWidth',2)
%     plot(t_,E0 + ft([dE,tau],t_),'--b','LineWidth',2)
%     plot(t_,true_params(1) + ft([true_params(2),true_params(3)],t_),':g','LineWidth',2)
% end
% figure;
% for i=1:3
%     subplot(1,3,i);
%     histogram(bias(:,i));
% end
