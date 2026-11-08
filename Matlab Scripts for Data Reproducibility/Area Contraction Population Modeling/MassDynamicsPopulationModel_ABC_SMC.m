%% SMC Sampler for Mass Dynamics Population
%{
    Last Updated: MJB 11/7/2026

    Calculates the posterior distribution of population model parameters
    for the mass dynamics model using Approximate Bayesian Computation 
    (synthetic likelihoods) Sequential Monte Carlo methods.
%}

clear all; close all; clc

addpath(cd+"\Sequential Monte Carlo Library")
addpath(cd+"\Population Models")

% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool("local", 8);  % Adjust the number of workers if needed
end

timestamp = datestr(now, 'yyyymmdd_HH00');

% Parameters related to saving the figures and data
to_save = 1;
if to_save
    output_fig_folder = "Mass Results" +string(timestamp)+ "\";
    if isfolder(output_fig_folder)
        % Remove the folder and all of its contents status safely
        rmdir(output_fig_folder, 's'); 
    end
    mkdir(output_fig_folder)
end

% setting for reproducibility
rng(1000000); 


%% Loading the Experimental data
root = "../Experimental Data/";

% loading the inventory table with all of the summary data
inventory_file = "Inventory Final.xlsx";
inventory = readtable(root + inventory_file,"ReadVariableNames",true);

% pulling out relevant columns of the inventory table (CrosslinkTime,
% CrosslinkTemp, Starting Mass, Final Mass)
data = inventory(:, [4, 3, 5, 6]);

rows27 = find(inventory.CrosslinkTemp == 27)';
rows37 = find(inventory.CrosslinkTemp == 37)';
rows0 = find(inventory.CrosslinkTemp == 0)'; % for samples that were tested immediately with no additional cross-linking

rows = [rows0, rows27, rows37];
data = data{rows,:}; % 27C and 37C Initial and Final mass data
data(:,2) = data(:,2)*60; % converting to minutes for agreement with the area data units

data(data(:,1)==0,1) = -1; 
data(data(:,1)>0,1) = data(data(:,1)>0,1)==37;
data_27 = data((data(:,1)==0)|(data(:,1)==-1),:);
data_37 = data((data(:,1)==1)|(data(:,1)==-1),:);

%% Setting up the priors

% Mean Starting Mass [g]
% centered on the value specified in the CAD
mu_M0_prior = makedist("Normal","mu",1, "sigma",1); 

% Mean Change in Mass [g] (assuming no average change)
mu_dM_prior = makedist("Normal","mu",0, "sigma",0.5);

% Standard Deviation of Mass Variables
% Half normal to enforce a strictly positive standard deviation. 0.5 so broad relative to mean
sig_M_prior = makedist("HalfNormal",'sigma',0.5);

% Correlation coefficient: fully uninformative
rho_M_prior = makedist("Uniform", "lower", -1, "upper", 1);


% Mean and Standard Deviation of Time Constant [min]
mean_tau_prior = makedist("Gamma", "a", 3, "b", 40);  % mean = 120 min, std ≈ 69 min
sig_tau_prior = makedist("Gamma", "a", 2, "b", 30);  % mean = 60 min, std ≈ 42 min

% Prior for the effect sizes. Truncated Normal to ensure that temperature
% doesn't change the sign of the CSA change
effect_min = -0.9999;
effect_max = 5;
effect_prior = truncate(makedist("Normal",0,1),effect_min,effect_max);


%% Fitting the Posterior Parameter Distribution

% geometric interpolation coefficients
gammas_ = linspace(0.0,1,10)'; 

% names of the parameters being estimated
parameter_labels = {"\mu_{M_0}","\mu_{\Delta M}",...
                    "\sigma_{M_0}","\sigma_{\Delta M}","\rho",...
                    "\mu_\tau","\sigma_\tau",...
                    "\alpha_\mu","\alpha_\sigma",...
                    "\beta_\mu","\beta_\sigma"};

% parameter names simplified for names
text_names = {"mu_M0","mu_dM",...
            "sigma_M0","sigma_dM","rho",...
            "mu_tau","sigma_tau",...
            "alpha_mu","alpha_sigma",...
            "beta_mu","beta_sigma"};

% population model class specialized for the Mass measurement experiments
model = Mass_PopulationModel(1);

mappings = {};  % No mappings because we have independent priors

priors = {mu_M0_prior,mu_dM_prior,...
          sig_M_prior, sig_M_prior, rho_M_prior,...
          mean_tau_prior, sig_tau_prior,...
          effect_prior,effect_prior,effect_prior,effect_prior};

% transforming parameters to bring bounded parameters to unbounded domain
transforms = {{3, LogTransform()};                              % standard deviation of initial CSA
              {4, LogTransform()};                              % standard deviation of change in CSA
              {5, LogitTransform(-1,1)};                        % correlation coefficient
              {6, LogTransform()};                              % mean time constant
              {7, LogTransform()};                              % std time constant
              {8, LogitTransform(effect_min,effect_max)};       % effects
              {9, LogitTransform(effect_min,effect_max)};
              {10, LogitTransform(effect_min,effect_max)};
              {11, LogitTransform(effect_min,effect_max)};};

% setting up the SMC sampler
sampler = SMC_Sampler(5000, model, priors,2.,data,gammas_,parameter_labels, mappings, transforms,0);

% generating the posterior
sampler.GeneratePosterior();
if to_save
    save(output_fig_folder+"SMC_sampler_results.mat","sampler")
end

%% Plotting the convergence and the summary of the posterior
figs = sampler.PlotConvergence();
for i = 1:length(figs)
    figure(figs{i});
    set(gca,'FontName','Calibri','FontSize',18)
    ylabel("")
    xlabel("")
    if to_save
        saveas(gcf, output_fig_folder + text_names{i} + "_Convergence.svg", "svg");
    end
end

[summary,lines] = sampler.PosteriorSummary();
if to_save
    writelines(lines, output_fig_folder + "PosteriorSummary.txt");
end

%% Comparison of the prior and the posterior

for i=1:length(priors)
    fig = figure('position',[100,100,400,200],'color','w');hold all
    ax = gca;
    ax.PositionConstraint = 'innerposition'; % Keeps the inner plot area fixed

    fig = sampler.PlotPrior(fig,i);         % plot the prior
    fig = sampler.PlotKDEstimate(fig,i);    % plot the kernel density estimate of the posterior

    % for the effects, add a verticle line at 0
    if i>length(priors)-4
        xline(0,'--k','LineWidth',2)
    end

    set(gca,'FontName','Calibri','FontSize',18)
    ylabel("")
    xlabel("")
    if to_save
        saveas(gcf, output_fig_folder + text_names{i} + ".svg", "svg");
    end
end

%% Median Posterior Simulations
theta_post = summary(:,1);      % getting the posterior median particle
model.setParams(theta_post);    % setting the model class to have the posterior median parameters

t = 0:5:360;       % time points at which to calculate the predicted trajectories
N_traj = 100000;    % number of trajectories to sample

% trajectories for the 27C samples
[trajectories,~] = model.rnd_trajectories(t,N_traj,0);
mean_27 = mean(trajectories,1);                 % mean
std_27 = std(trajectories,[],1);                % standard deviation
lower_CI_27 = quantile(trajectories,0.025,1);   % upper and lower bounds of the 95% credible interval
upper_CI_27 = quantile(trajectories,0.975,1);

% trajectories for the 37C samples
[trajectories,particles] = model.rnd_trajectories(t,N_traj,1);
mean_37 = mean(trajectories,1);
std_37 = std(trajectories,[],1);
lower_CI_37 = quantile(trajectories,0.025,1);
upper_CI_37 = quantile(trajectories,0.975,1);

%% Posterior Predictive Check

color_27 = [0,192,244]/255;
color_37 = [235,72,0]/255; 

% 27C Data
ind_27 = find((data(:,1)==0)|(data(:,1)==-1));
f1 = figure("Position",[100,100,500,300],'Color','w'); hold all
fill([t,flip(t)],[upper_CI_27,flip(lower_CI_27)],color_27,'EdgeColor','none','FaceAlpha',0.2);
fill([t,flip(t)],[mean_27-std_27,flip(mean_27 + std_27)],color_27,'EdgeColor','none','FaceAlpha',0.5);

plot(data(ind_27,2),data(ind_27,4),'ok','LineWidth',1)

plot(t,mean_27,'-','color',0.6*color_27,'LineWidth',3)
xlabel("Time [min]")
ylabel("Mass [g]")
% ylim([30,415])
set(gca,'FontName','Calibri','FontSize',14)
if to_save
    saveas(gcf, output_fig_folder + "GroupPPC_27C.svg", "svg");
end


% 37C Data
ind_37 = find((data(:,1)==1)|(data(:,1)==-1));
f1 = figure("Position",[100,100,500,300],'Color','w'); hold all
fill([t,flip(t)],[upper_CI_37,flip(lower_CI_37)],color_37,'EdgeColor','none','FaceAlpha',0.2);
fill([t,flip(t)],[mean_37-std_37,flip(mean_37 + std_37)],color_37,'EdgeColor','none','FaceAlpha',0.5);

plot(data(ind_37,2),data(ind_37,4),'ok','LineWidth',1)

plot(t,mean_37,'-','color',0.6*color_37,'LineWidth',3)
xlabel("Time [min]")
ylabel("Mass [g]")
% ylim([30,415])
set(gca,'FontName','Calibri','FontSize',14)
if to_save
    saveas(gcf, output_fig_folder + "GroupPPC_37C.svg", "svg");
end


% Comparison of the groups
f1 = figure("Position",[100,100,500,300],'Color','w'); hold all
fill([t,flip(t)],[mean_27-std_27,flip(mean_27 + std_27)],color_27,'EdgeColor','none','FaceAlpha',0.5);
plot(t,mean_27,'-','color',0.6*color_27,'LineWidth',3)
fill([t,flip(t)],[mean_37-std_37,flip(mean_37 + std_37)],color_37,'EdgeColor','none','FaceAlpha',0.5);
plot(t,mean_37,'-','color',0.6*color_37,'LineWidth',3)
xlabel("Time [min]")
ylabel("Mass [g]")
% ylim([30,400])
set(gca,'FontName','Calibri','FontSize',14)
if to_save
    saveas(gcf, output_fig_folder + "GroupComparison.svg", "svg");
end

%% Plotting the Joint Posteriors

% M, dM:
fig1 = figure("Position",[100,100,750,600],'Color','w');
sampler.PlotJointPosterior(fig1,[1,2,3,4,5])
if to_save
    saveas(fig1, output_fig_folder + "M_dM_JointPrior.svg", "svg");
end

% dM, effects:
fig2 = figure("Position",[100,100,750,600]*0.8,'Color','w');
sampler.PlotJointPosterior(fig2,[2,4,8,9])
if to_save
    saveas(fig2, output_fig_folder + "dM_Effects_JointPrior.svg", "svg");
end

% tau, effects:
fig3 = figure("Position",[100,100,750,600]*0.8,'Color','w');
sampler.PlotJointPosterior(fig3,[6,7,10,11])
if to_save
    saveas(fig3, output_fig_folder + "tau_Effects_JointPrior.svg", "svg");
end
