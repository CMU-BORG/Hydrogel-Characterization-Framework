%% Script to regenerate plots for previous results
%{
    Last Updated: MJB 11/24/2026
    
    This script allows for the results of the SMC samplers to be generated
    without having to rerun the simulations.
%}

clear all; close all; clc

cur_dir = split(cd,"\");
root_dir = join(cur_dir(1:end-2),"\");

addpath(root_dir+"\Sequential Monte Carlo Library")
addpath(root_dir+"\Population Models")

sampler_file_name = "SMC_sampler_results.mat";

timestamp = datestr(now, 'yyyymmdd_HH00');

%% Area Results
area_folder = "Area Results20261108_0000\";
area_data = load(area_folder + sampler_file_name).sampler;

text_names = {"mu_A0","mu_dA",...
            "sigma_A0","sigma_dA","rho",...
            "mu_tau","sigma_tau",...
            "alpha_mu","alpha_sigma",...
            "beta_mu","beta_sigma"};

PlotResults(area_data,text_names,area_folder,"Cross-Sectional Area [mm^2]");

writelines(sprintf("~~~~~ NOTE ~~~~~\nFigures regenerated on %s.",timestamp),area_folder+"README.txt")

close all
%% Mass Results
mass_folder = "Mass Results20261108_0000\";
mass_data = load(mass_folder + sampler_file_name).sampler;

text_names = {"mu_M0","mu_dM",...
            "sigma_M0","sigma_dM","rho",...
            "mu_tau","sigma_tau",...
            "alpha_mu","alpha_sigma",...
            "beta_mu","beta_sigma"};

PlotResults(mass_data,text_names,mass_folder,"Mass [g]");

writelines(sprintf("~~~~~ NOTE ~~~~~\nFigures regenerated on %s.",timestamp),mass_folder+"README.txt")

close all
%% Stiffness Results
stiff_folder = "Stiffness Results20261108_0100\";
stiff_data = load(stiff_folder + sampler_file_name).sampler;

text_names = {"mu_E0","mu_dE",...
            "sigma_E0","sigma_dE","rho",...
            "mu_tau","sigma_tau",...
            "alpha_mu","alpha_sigma",...
            "beta_mu","beta_sigma"};

PlotResults(stiff_data,text_names,stiff_folder,"Effective Elastic Modulus [kPa]");

writelines(sprintf("~~~~~ NOTE ~~~~~\nFigures regenerated on %s.",timestamp),stiff_folder+"README.txt")

close all



function PlotResults(sampler, text_names, output_fig_folder, variable_label,to_save)
arguments
    sampler 
    text_names 
    output_fig_folder 
    variable_label 
    to_save = 1;
end
%{
    This function regenerates all of the plots that we need from the
    original run of the SMC sampler for a given dataset. No additional
    simulation are run, and all results are equivalent to previously
    generated results (save for changes in the plotting style).
%}

% ~~~~~~~~ Remove previous SVGs from the folder ~~~~~~~~ 
delete(output_fig_folder + "*.svg")

% ~~~~~~~~ Plotting the convergence and the summary of the posterior ~~~~~~~~ 
figs = sampler.PlotConvergence();
for i = 1:length(figs)
    figure(figs{i});
    set(gca,'FontName','Calibri','FontSize',18)
    ylabel("")
    xlabel("")
    if to_save
        set(gcf, 'Renderer', 'painters'); 
        saveas(gcf, output_fig_folder + text_names{i} + "_Convergence.svg", "svg");
    end
end


% ~~~~~~~~ Comparison of the prior and the posterior ~~~~~~~~ 
priors = sampler.priors;
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
        set(gcf, 'Renderer', 'painters'); 
        saveas(gcf, output_fig_folder + text_names{i} + ".svg", "svg");
    end
end

summary = sampler.PosteriorSummary();
model = sampler.model;

% ~~~~~~~~ Median Posterior Simulations ~~~~~~~~ 
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
[trajectories,~] = model.rnd_trajectories(t,N_traj,1);
mean_37 = mean(trajectories,1);
std_37 = std(trajectories,[],1);
lower_CI_37 = quantile(trajectories,0.025,1);
upper_CI_37 = quantile(trajectories,0.975,1);

% ~~~~~~~~ Posterior Predictive Check ~~~~~~~~ 
data = sampler.data;

color_27 = [0,192,244]/255;
color_37 = [235,72,0]/255; 

% 27C Data
ind_27 = find((data(:,1)==0)|(data(:,1)==-1));
f1 = figure("Position",[100,100,500,300],'Color','w'); hold all
fill([t,flip(t)],[upper_CI_27,flip(lower_CI_27)],color_27,'EdgeColor','none','FaceAlpha',0.2);
fill([t,flip(t)],[mean_27-std_27,flip(mean_27 + std_27)],color_27,'EdgeColor','none','FaceAlpha',0.5);

if strcmp(class(sampler.model),'PopulationModel') % need this specific check instead of 'isa' because 'isa' returns true for parent class as well
    % data for area data
    for i=1:length(ind_27)
        theta = data(ind_27(i),2:end);
        y_ = (theta(1) + theta(2)) - theta(2)*exp(-t/theta(3));
        plot(t,y_,'--k','LineWidth',0.1)
    end
else
    % data for mass and stiffness data
    plot(data(ind_27,2),data(ind_27,end),'ok','LineWidth',1)
end

plot(t,mean_27,'-','color',0.6*color_27,'LineWidth',3)
xlabel("Time [min]")
ylabel(variable_label)
ylim([0.85*min([min(lower_CI_27), min(lower_CI_37)]), 1.15*max([max(upper_CI_27), max(upper_CI_37)])])
set(gca,'FontName','Calibri','FontSize',14)
if to_save
    set(gcf, 'Renderer', 'painters'); 
    saveas(gcf, output_fig_folder + "GroupPPC_27C.svg", "svg");
end


% 37C Data
ind_37 = find((data(:,1)==1)|(data(:,1)==-1));
f1 = figure("Position",[100,100,500,300],'Color','w'); hold all
fill([t,flip(t)],[upper_CI_37,flip(lower_CI_37)],color_37,'EdgeColor','none','FaceAlpha',0.2);
fill([t,flip(t)],[mean_37-std_37,flip(mean_37 + std_37)],color_37,'EdgeColor','none','FaceAlpha',0.5);

if strcmp(class(sampler.model),'PopulationModel')
    % data for area data
    for i=1:length(ind_37)
        theta = data(ind_37(i),2:end);
        y_ = (theta(1) + theta(2)) - theta(2)*exp(-t/theta(3));
        plot(t,y_,'--k','LineWidth',0.1)
    end
else
    % data for mass and stiffness data
    plot(data(ind_37,2),data(ind_37,end),'ok','LineWidth',1)
end

plot(t,mean_37,'-','color',0.6*color_37,'LineWidth',3)
xlabel("Time [min]")
ylabel(variable_label)
ylim([0.85*min([min(lower_CI_27), min(lower_CI_37)]), 1.15*max([max(upper_CI_27), max(upper_CI_37)])])
set(gca,'FontName','Calibri','FontSize',14)
if to_save
    set(gcf, 'Renderer', 'painters'); 
    saveas(gcf, output_fig_folder + "GroupPPC_37C.svg", "svg");
end

% Comparison of the groups
f1 = figure("Position",[100,100,500,300],'Color','w'); hold all
fill([t,flip(t)],[mean_27-std_27,flip(mean_27 + std_27)],color_27,'EdgeColor','none','FaceAlpha',0.5);
plot(t,mean_27,'-','color',0.6*color_27,'LineWidth',3)
fill([t,flip(t)],[mean_37-std_37,flip(mean_37 + std_37)],color_37,'EdgeColor','none','FaceAlpha',0.5);
plot(t,mean_37,'-','color',0.6*color_37,'LineWidth',3)
xlabel("Time [min]")
ylabel(variable_label)
ylim([0.9*min([min(mean_27-std_27),min(mean_37-std_37)]),1.1*max([max(mean_27+std_27),max(mean_37+std_37)])])
set(gca,'FontName','Calibri','FontSize',14)
if to_save
    set(gcf, 'Renderer', 'painters'); 
    saveas(gcf, output_fig_folder + "GroupComparison.svg", "svg");
end


% ~~~~~~~~ Plotting the Joint Posteriors ~~~~~~~~ 

% Control Model Parameters:
fig1 = figure("Position",[100,100,750,600],'Color','w');
sampler.PlotJointPosterior(fig1,[1,2,3,4,5]);
if to_save
    set(fig1, 'Renderer', 'painters'); 
    saveas(fig1, output_fig_folder + "ControlModel_JointPrior.svg", "svg");
end

% Change in Variable, effects:
fig2 = figure("Position",[100,100,750,600]*0.8,'Color','w');
sampler.PlotJointPosterior(fig2,[2,4,8,9]);
if to_save
    set(fig2, 'Renderer', 'painters'); 
    saveas(fig2, output_fig_folder + "Change_Effects_JointPrior.svg", "svg");
end

% tau, effects:
fig3 = figure("Position",[100,100,750,600]*0.8,'Color','w');
sampler.PlotJointPosterior(fig3,[6,7,10,11]);
if to_save
    set(fig3, 'Renderer', 'painters'); 
    saveas(fig3, output_fig_folder + "Tau_Effects_JointPrior.svg", "svg");
end

end