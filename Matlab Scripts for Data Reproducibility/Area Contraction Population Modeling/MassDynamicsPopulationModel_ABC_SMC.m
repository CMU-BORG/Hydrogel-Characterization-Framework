%% SMC Sampler for Mass Dynamics Population
%{


%}

clear all; close all; clc

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
% CrosslinkTemp, Scale_pxpmm)
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
%%
figure; hold all
% plot(0*data(:,2),data(:,3),'^k')
% plot(data_27(:,2),data_27(:,4),'sr');
% plot(data_37(:,2),data_37(:,4),'ob');
% errorbar(0*data(:,2),median(data(:,3)), quantile(data(:,3),[0.25]) - median(data(:,3)), quantile(data(:,3),[0.75]) - median(data(:,3)),'sr',"LineWidth",2)

plot(0*data(:,2),0*data(:,3),'sk')

for i=60*[2,4,6]
    inds = find(data_27(:,2)==i);
    data_i = data_27(inds,4)-data_27(inds,3);
    errorbar(data_27(inds,2),median(data_i), quantile(data_i,[0.25]) - median(data_i), quantile(data_i,[0.75]) - median(data_i),'sr',"LineWidth",2)
    
    inds = find(data_37(:,2)==i);
    data_i = data_37(inds,4)-data_37(inds,3);
    errorbar(data_37(inds,2),median(data_i), quantile(data_i,[0.25]) - median(data_i), quantile(data_i,[0.75]) - median(data_i),'sb',"LineWidth",2)
    
end
%% Setting up the priors

mu_M0_mean = 1; % [g] (centered on the value specified in the CAD)
mu_M0_std = 1;  % [1] (variance is known from printer calibration to not be too large)
mu_M0_prior = makedist("Normal","mu",mu_M0_mean, "sigma",mu_M0_std);


mu_dM_mean = 0; % [g] assume no change in mass
mu_dM_std = 1;  % [g] 
mu_dM_prior = makedist("Normal","mu",mu_dM_mean, "sigma",mu_dM_std);

sig_M_lower = 0.00001;      % [g] lowest possible value of sigma  
sig_M_upper = 1;    % [g] highest possible values of sigma
sig_M_prior = makedist("Uniform","lower",sig_M_lower,"upper",sig_M_upper);

% CORRELATION COEFFICIENT PRIOR: 
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
effect_prior = makedist("Uniform","lower",-0.9999,"upper",5); 

priors = {mu_M0_prior,mu_dM_prior,...
          sig_M_prior, sig_M_prior, rho_A_prior,...
          mean_tau_prior, sig_tau_prior,...
          effect_prior,effect_prior,effect_prior,effect_prior};

parameter_labels = {"\mu_{M_0}","\mu_{\Delta M}",...
                    "\sigma_{M_0}","\sigma_{\Delta M}","\rho",...
                    "\mu_\tau","\sigma_\tau",...
                    "\alpha_\mu","\alpha_\sigma",...
                    "\beta_\mu","\beta_\sigma"};


%% Step 1: Estimate the initial mass distribution using standard SMC
gammas_ = linspace(0,1,10)';

% Getting an M0 distribution using all of the A0 data
model_norm_prior = GaussianModel(0);
M0_prior = SMC_Sampler(5000,model_norm_prior,priors([1,3]),2.38,data(:,3),gammas_,{"\mu_{M_0}","\sigma_{M_0}"},0,{},1);
M0_prior.GeneratePosterior();
M0_prior.PlotConvergence();

%% Step 2: Estimating the Population model for 27C and 37C separately

mass_pop_model_one_temp = Mass_PopulationModel(0);

N_gens = 20;
gammas = [linspace(1,0.,N_gens)]';
gamma_min = 0.1;
gammas = (1-gamma_min)*(3*gammas.^2 - 2*gammas.^3) + gamma_min;
% gammas = zeros(N_gens,1,4);
% for i=1:4
%     gammas(:,1,i) = 3*gammas_mu.^2 - 2*gammas_mu.^3;
% end

% 27 C data
time_vals = unique(data(:,2));
data_cell = cell(length(time_vals),1);
for i=1:length(time_vals)
    data_i = data_37(data_37(:,2)==time_vals(i),:);
    data_cell{i} = {length(data_i), [time_vals(i), 0], data_i(:,[3,4])};
end

% constructing priors and mappings
priors_ = {M0_prior, priors{[2,4,5,6,7]} };
mappings_ = {1,[1,3];   
             2,[2];
             3,[4];
             4,[5];
             5,[6];
             6,[7]};

data_27_prior = SMC_Sampler(5000, mass_pop_model_one_temp, priors_, 2.38, data_cell, gammas, parameter_labels(1:7), 1, mappings_, 1);
% data_27_prior.S_sims = 50;
data_27_prior.GeneratePosterior();
data_27_prior.PlotConvergence();