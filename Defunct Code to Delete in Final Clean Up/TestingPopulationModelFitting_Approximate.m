%% Testing the Synthetic Likelihood Approach to Fitting the Population Model
%{

%}

clear all; close all; clc

addpath(cd+"\Sequential Monte Carlo Library")
addpath(cd+"\Population Models")


% Start parallel pool if not already runninger
if isempty(gcp('nocreate'))
    parpool("local", 8);  % Adjust the number of workers if needed
end

timestamp = datestr(now, 'yyyymmdd_HHMM');

%% Generating Synthetic Data using the Population Model

N_samps_per_group = 50; % use a lot of data to make sure its right
t_groups = 0:30:360;

data = zeros(N_samps_per_group*length(t_groups), 3);

model_true = Stiffness_PopulationModel(0);

% rng(1000000); % setting for reproducibility across runs
mu_E0_true = 150; sig_E0_true = 20;
mu_dE_true = 50;  sig_dE_true = 15;
rho_true = -0.;
mu_tau_true = 100; sig_tau_true = 20;
theta_true = [mu_E0_true,mu_dE_true,sig_E0_true,sig_dE_true,rho_true,mu_tau_true,sig_tau_true];

model_true.setParams(theta_true);

data_groups = cell(size(t_groups));
for i=1:length(t_groups)
    inds = 1+N_samps_per_group*(i-1):N_samps_per_group*i;
    data(inds,1) = -1*(t_groups(i)==0);
    data(inds,2) = t_groups(i);
    data(inds,3) = model_true.rnd_trajectories(t_groups(i),N_samps_per_group,0);
    data_groups{i} = data(inds,3);
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

sampler = SMC_Sampler(2000, model, priors,0.5,data,gammas_,parameter_labels, mappings_, transforms,0);
% 

%% Histograms
sampler.GeneratePosterior();

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
% theta_post = summary(:,1);

t_ = 0:30:360;
N_traj = 100000;

% model.setParams(theta_post);
% [trajectories,particles] = model.rnd_trajectories(t_,N_traj,0);
model.setParams(theta_true);
[trajectories_true,particles] = model.rnd_trajectories(t_,N_traj,0);

figure("Position",[100,100,600,400],"Color","w"); hold all

plot(data(:,2),data(:,3),'ok')
for i=1:length(t_groups)
    inds = 1+N_samps_per_group*(i-1):N_samps_per_group*i;
    data_ = data(inds,3);
    errorbar(t_groups(i),median(data_), -std(data_), std(data_),'sg',"LineWidth",2)
end

% mean_27 = mean(trajectories,1);
% std_27 = std(trajectories,[],1);
% errorbar(t_,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','r','Color','r')

mean_27 = mean(trajectories_true,1);
std_27 = std(trajectories_true,[],1);
errorbar(t_,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','b','Color','b')

% 
% %% Testing the Assumption of Mean Dynamics
% 
% ft = @(x,t_) x(1,:).*(1 - exp(-t_./x(2,:))); 
% % t_ = linspace(0,400,100);
% 
% data_ = reshape(data(:,3),[6,4]);
% 
% % model.setParams(sampler.SamplePrior(1));
% model.setParams(theta_true); %.*(1. + 0.1*(rand(size(theta_true))-0.5)));
% 
% replicates = 2000;
% RMSE_mean_true = zeros(replicates,1);
% RMSE_mean_data = zeros(replicates,1);
% RMSE_std_true = zeros(replicates,1);
% RMSE_std_data = zeros(replicates,1);
% for j=1:replicates
%     t = 0:120:360;
%     N_traj = 6;
% 
%     particles = zeros(N_traj*length(t),3);
%     trajectories_true = zeros(N_traj,length(t));
% 
%     for i=1:length(t)
%         [trajectories_true(:,i),particles(1+(i-1)*N_traj:i*N_traj,:)] = model.rnd_trajectories(t(i),N_traj,0);
%     end
% 
%     % what is we bootstrap the relaxation
%     % mean_trajectory = mean(trajectories_true,1);
%     N_boot = 100;
%     % E0 = mean_trajectory(1);
%     % [dE, tau] = model.FitRelaxation(t,mean_trajectory' - E0);
% 
%     inds = randi(N_traj,length(t),N_boot);
%     E_est_data = 0*inds;
%     E_est_sim = 0*inds;
%     for m=1:length(t)
%         E_est_sim(m,:) = trajectories_true(inds(m,:),m);
%         E_est_data(m,:) = data_(inds(m,:),m);
%     end
%     [dE_sim, tau_sim] = model.FitRelaxation(t,(E_est_sim - E_est_sim(1,:)));
%     [dE_data, tau_data] = model.FitRelaxation(t,(E_est_data - E_est_data(1,:)));
%     % plot(t_,E_est(1) + ft([dE,tau],t_),'--b','LineWidth',0.01);
%     E_boot_sim = E_est_sim(1,:) + ft([dE_sim;tau_sim],t_');
%     E_boot_data = E_est_data(1,:) + ft([dE_data;tau_data],t_');
% 
% 
%     % fit_params = [E0,dE,tau];
%     % true_params = [mean(particles(:,1)),mean(particles(:,2)),1/mean(1./particles(:,3))];
%     % bias_j = 100*(true_params - fit_params)./true_params;
%     % bias(j,:) = bias_j;
% 
%     % plot(t,mean_trajectory,'-r','LineWidth',2)
% 
%     mean_total_var_true = var(E_boot_sim,[],2) + std_27'.^2;
%     mean_total_var_data = var(E_boot_sim,[],2) + var(E_boot_data,[],2);
% 
%     % plot(t_,true_params(1) + ft([true_params(2),true_params(3)],t_),':g','LineWidth',2)
%     RMSE_mean_true(j) = (mean( ((mean(E_boot_sim,2) - mean_27')).^2./mean_total_var_true )) ;
%     RMSE_mean_data(j) = (mean( ((mean(E_boot_sim,2) - mean(E_boot_data,2))).^2./mean_total_var_data ));
% 
%     RMSE_std_true(j) = (mean( ((std(E_boot_sim,[],2) - std_27')).^2 ));
%     RMSE_std_data(j) = (mean( ((std(E_boot_sim,[],2) - std(E_boot_data,[],2))).^2 ));
% end
% 
% figure; hold all
% plot(data(:,2),data(:,3),'ok')
% errorbar(t_,mean(E_boot_data,2),std(E_boot_data,[],2),std(E_boot_data,[],2),'vertical','Marker','^','MarkerFaceColor','b','Color','b')
% errorbar(t_,mean(E_boot_sim,2),std(E_boot_sim,[],2),std(E_boot_sim,[],2),'vertical','Marker','o','MarkerFaceColor','g','Color','g')
% errorbar(t_,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','r','Color','r')
% 
% ksi2_pdf = makedist("Gamma","a",length(t_)/2,"b",2);
% 
% figure; 
% subplot(1,2,1);hold all
% histogram(RMSE_mean_data,'Normalization','pdf','DisplayName',"Data")
% histogram(RMSE_mean_true,'Normalization','pdf','DisplayName',"True")
% ksi2_pdf.plot()
% legend()
% subplot(1,2,2);hold all
% histogram(RMSE_std_data,'Normalization','pdf','DisplayName',"Data")
% histogram(RMSE_std_true,'Normalization','pdf','DisplayName',"True")
% legend()
% 
% fprintf("Average Mean RMSE Compared to Data: %.3f +/- %.3f\n",mean(RMSE_mean_data),std(RMSE_mean_data))
% fprintf("Average Mean RMSE Compared to True: %.3f +/- %.3f\n",mean(RMSE_mean_true),std(RMSE_mean_true))
% 
% fprintf("Average STD RMSE Compared to Data: %.3f +/- %.3f\n",mean(RMSE_std_data),std(RMSE_std_data))
% fprintf("Average STD RMSE Compared to True: %.3f +/- %.3f\n",mean(RMSE_std_true),std(RMSE_std_true))


%% Seeing how we can look at the distribution of the errors
t_ = 0:30:360;
ft = @(x,t_) x(1,:).*(1 - exp(-t_./x(2,:))); 
% t_ = linspace(0,400,100);

data_ = reshape(data(:,3),[6,4]);

% model.setParams(sampler.SamplePrior(1));
model.setParams(theta_true); %.*(1. + 0.1*(rand(size(theta_true))-0.5)));

N_boot = 500; N_traj = 6;
inds = randi(N_traj,length(t_groups),N_boot);
    
% data bootstrap
E_est_data = 0*inds;
for m=1:length(t_groups)
    E_est_data(m,:) = data_(inds(m,:),m);
end
[dE_data, tau_data] = model.FitRelaxation(t_groups,(E_est_data - E_est_data(1,:)));
E_boot_data = 0*E_est_data(1,:) + ft([dE_data;tau_data],t_');

data_mean = mean(E_boot_data,2);
data_std = std(E_boot_data,[],2);

figure; hold all
plot(t_groups,data_,'ok')
plot(t_,E_boot_data,':k','LineWidth',0.001)
errorbar(t_,data_mean,data_std,data_std,'vertical','Marker','^','MarkerFaceColor','b','Color','b')
errorbar(t_,mean_27,std_27,std_27,'vertical','Marker','s','MarkerFaceColor','r','Color','r')


replicates = 200;
t = 0:120:360;
N_traj = 6;

mean_sim_data = zeros(replicates,length(t_));
std_sim_data = zeros(replicates,length(t_));

% for j=1:replicates
    % trajectories_true = zeros(N_traj,length(t));
    
    [trajectories_true,~] = model.rnd_trajectories(t,N_traj*replicates,0);
    trajectories_true = reshape(trajectories_true,[replicates,N_traj,length(t)]);
    
    % what is we bootstrap the relaxation
    E_est_sim = zeros(replicates,length(t_groups),N_boot);
    for m=1:length(t_groups)
        E_est_sim(:,m,:) = trajectories_true(:,inds(m,:),m);
    end
    E_est_sim = permute(E_est_sim,[2,3,1]);
    E_all = E_est_sim;
    E_est_sim = reshape(E_est_sim,size(E_est_sim,1),[]);

    [dE_sim, tau_sim] = model.FitRelaxation(t_groups,(E_est_sim - E_est_sim(1,:)));
    E_boot_sim = 0*E_est_sim(1,:) + ft([dE_sim;tau_sim],t_');
    E_boot_sim = reshape(E_boot_sim,[length(t_),N_boot,replicates]);
    mean_sim_data = permute(mean(E_boot_sim,2),[1,3,2]);
    std_sim_data = permute(std(E_boot_sim,[],2),[1,3,2]);
% end


figure('Position',[100,100,600,300]); 
subplot(1,2,1); hold all
plot(t_,mean_sim_data,':k','LineWidth',0.001)
errorbar(t_,mean(mean_sim_data,2),std(mean_sim_data,[],2),std(mean_sim_data,[],2),'vertical','Marker','^','MarkerFaceColor','k','Color','k','LineWidth',2)
plot(t_,data_mean,'-b','LineWidth',2)

subplot(1,2,2); hold all
plot(t_,std_sim_data,':k','LineWidth',0.001)
errorbar(t_,mean(std_sim_data,2),std(std_sim_data,[],2),std(std_sim_data,[],2),'vertical','Marker','^','MarkerFaceColor','k','Color','k','LineWidth',2)
plot(t_,data_std,'-b','LineWidth',2)

p_vec = normpdf(data_mean(2:end),mean(mean_sim_data(2:end,:),2),std(mean_sim_data(2:end,:),[],2));
p_ = exp(mean(log(p_vec)));
fprintf("Probability of mean vector: %.3f\n",p_)

p_vec = normpdf(data_std(2:end),mean(std_sim_data(2:end,:),2),std(std_sim_data(2:end,:),[],2));
p_ = exp(mean(log(p_vec)));
fprintf("Probability of std vector: %.3f\n",p_)