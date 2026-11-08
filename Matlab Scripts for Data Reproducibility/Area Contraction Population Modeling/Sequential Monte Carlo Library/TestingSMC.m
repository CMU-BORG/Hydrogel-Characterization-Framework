%% Testing SMC Sampler 
%{
    Last Updated: MJB 11/7/2026

    This script validates the algorithm implementation of the Sequential
    Monte Carlo method using known solutions or integrable posteriors.
%}

clear all; close all; clc;

% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool("local", 8);  % Adjust the number of workers if needed
end

%% Testing Using Bernoulli Model and Conjugate Prior
%{
    As a first step, we can estimate the probability parameter of a
    Bernoulli model with a Beta prior. A Beta distribution is a conjugate
    prior, so we can calculate the exact analytical solution and use it to
    validate the SMC estimate.
%}

bernoulli_model = BernoulliModel();    % generating model

N_data = 20;    % number of datapoints to use in the dataset
p_true = 0.75;  % true value of the model parameter for generating the data

bernoulli_model.setParams([p_true]);        % set the true parameter and generate the calibration data
X_data_true = bernoulli_model.rnd(N_data);
figure; histogram(X_data_true,'Normalization','probability','BinWidth',0.333);  % plot the calibration data

% Set up the conjugate prior and the resulting posterior
a_prior = 1; b_prior = 1;                               % hyper parameters of the prior
priors = {makedist("Beta","a",a_prior,"b",b_prior)};    % using a Beta prior so we can test the SMC algorithm against a known conjugate prior
posterior = makedist("Beta","a",a_prior+sum(X_data_true),"b",b_prior+N_data-sum(X_data_true));    % true posterior based on conjugacy

gammas = min(linspace(0,1.,10),1)';         % gamma values to use for the SMC algorithm
transforms = {{1, LogitTransform(0,1)};};   % applying a Logit transformation to the parameter (bound between 0 and 1) to bring it into the range (-inf, inf)

sampler = SMC_Sampler(1000,bernoulli_model,priors,2.38,X_data_true,gammas,["p"],{},transforms); % set up the SMC sampler and run the calibration
sampler.GeneratePosterior();

% Plotting the partial posterior across the generations and comparing to
% the true posterior distribution
figs = sampler.PlotConvergence();
figure(figs{1}); hold all
yline(p_true,'--k','LineWidth',2)
yline(posterior.mean(),'-b','LineWidth',2)
yline(posterior.icdf([0.025,0.975]),':b','LineWidth',2)

% Comparing to the known PDFs
fig_post = sampler.PlotPosteriorHistograms();
x_ = linspace(0,1,300);
figure(fig_post{1}); hold all
plot(x_,priors{1}.pdf(x_),'--','Color','r','LineWidth',2)
plot(x_,posterior.pdf(x_),'-','Color','b','LineWidth',2);
plot(x_,sampler.pdf(x_'),':g','LineWidth',2)
histogram(sampler.random(5000),'FaceAlpha',0.2,'Normalization','pdf')

% Comparing to the kernel density estimates of the PDFs
fig_post_KDE = figure('Position',[100,100,400,200],'Color','w');
fig_post_KDE = sampler.PlotPrior(fig_post_KDE,1);
fig_post_KDE = sampler.PlotKDEstimate(fig_post_KDE,1);




%% Testing with Gamma for multi-variable fitting
%{
    Next, we can test fitting multiple variables using a Gamma generating
    model. There is not a simple conjugate prior that we can use, so we can
    numerically integrate the posterior (using traditional integration
    techniques) to validate the sampler.
%}

gamma_model = GammaModel(0);   % setting up the generating 

N_data = 30;    % number of datapoints to use in the dataset
alpha_true = 1;  % true value of the model parameters for generating the data
beta_true = 10;

gamma_model.setParams([alpha_true,beta_true]);  % setting the parameters and generating the data
X_data = gamma_model.rnd(N_data);
figure; histogram(X_data,'Normalization','probability');

% hyper parameters for a Uniform prior over both parameters
alpha_prior_a = 0; alpha_prior_b = 5; 
beta_prior_a = 0; beta_prior_b = 50; 
priors = {makedist("Uniform","lower",alpha_prior_a,"upper",alpha_prior_b),
          makedist("Uniform","lower",beta_prior_a,"upper",beta_prior_b)}; 


% numerically calculating the true posterior on a grid over the priors
N_points = 100;
alpha_test = linspace(alpha_prior_a,alpha_prior_b,N_points+1);
alpha_test(1) = []; % removing 0 to avoid issues with the Gamma distribution
beta_test = linspace(beta_prior_a,beta_prior_b,N_points+1);
beta_test(1) = [];

[A,B] = meshgrid(alpha_test,beta_test); % generating the integration grid
post_ = 0*A;    % preallocating
for i=1:N_points
    for j=1:N_points
        % loop over all of the points in the grid, set the parameters for
        % the model, and calculate the un-normalized posterior
        gamma_model.setParams([A(i,j),B(i,j)]);
        post_(i,j) = exp(gamma_model.log_likelihood(X_data) + log(priors{1}.pdf(A(i,j))) + log(priors{2}.pdf(B(i,j))) ); 
    end
end

% numerically integrate the posterior to normalize
I = trapz(beta_test, trapz(alpha_test, post_, 2));
post_ = post_ / I;  % numerically integrated posterior

% Using numerical integration to get the marginal posteriors for each
% parameter and the marginal CDFs
f_alpha = trapz(beta_test, post_, 1);
f_beta = trapz(alpha_test, post_, 2);
F_alpha = cumtrapz(alpha_test,f_alpha); % CDF
F_beta = cumtrapz(beta_test,f_beta); % CDF

% using interpolation on the CDFs to calculate the 95% credible intervals
% and median on the parameters from the numerical estimates
val_Fa = (F_alpha>0.001)&(F_alpha<0.999);
val_Fb = (F_beta>0.001)&(F_beta<0.999);
alpha_inv_CDF = @(p) interp1(F_alpha(val_Fa),alpha_test(val_Fa),p);
beta_inv_CDF = @(p) interp1(F_beta(val_Fb),beta_test(val_Fb),p);
alpha_CI = alpha_inv_CDF([0.025,0.975]);
beta_CI = beta_inv_CDF([0.025,0.975]);
alpha_med = alpha_inv_CDF(0.5);
beta_med = beta_inv_CDF(0.5);

% plotting the numerically integrated posteriors
figure; hold all
plot(alpha_test,f_alpha,'-r','LineWidth',2);
plot(beta_test,f_beta,'-b','LineWidth',2);

% Setting up the sampler 
gammas = linspace(0,1.,10)';   % gamma values to use for the SMC algorithm
transforms = { {1,LogTransform()}, {2,LogTransform()} };

sampler = SMC_Sampler(1000,gamma_model,priors,2.38,X_data,gammas,["\alpha","\beta"],{},transforms);
sampler.GeneratePosterior();

% getting the marginal pdfs to compare with the numerically integrated
% values
fa_post = sampler.marginal_pdf(1,alpha_test');
fb_post = sampler.marginal_pdf(2,beta_test');

% calculating the parameters we would expect from MLE
mu_data = mean(X_data);
sig_data = std(X_data);
alpha_data = (mu_data/sig_data)^2;
beta_data = sig_data^2 / mu_data;

% generating posterior samples to validate the resampler
post_resamps = sampler.random(5000);

% plotting the posterior histograms
fig_post = sampler.PlotPosteriorHistograms();

% plotting alpha results against the true values and the numerically
% integrated posteriors
x_ = linspace(0,50,300);
figure(fig_post{1}); hold all
xline(alpha_true,'--k','LineWidth',2);
xline(alpha_data,'--g','LineWidth',2);
plot(0.05*x_,priors{1}.pdf(0.05*x_),'--','Color','r','LineWidth',2)
plot(alpha_test,f_alpha,'--b','LineWidth',2)
plot(alpha_test,fa_post,':g','LineWidth',2)
histogram(post_resamps(:,1),'FaceAlpha',0.2,'Normalization','pdf')

% repeating for beta
figure(fig_post{2}); hold all
xline(beta_true,'--k','LineWidth',2);
xline(beta_data,'--g','LineWidth',2);
plot(x_,priors{2}.pdf(x_),'--','Color','r','LineWidth',2)
plot(beta_test,f_beta,'--b','LineWidth',2)
plot(beta_test,fb_post,':g','LineWidth',2)
histogram(post_resamps(:,2),'FaceAlpha',0.2,'Normalization','pdf')


% Comparing to the kernel density estimates of the PDFs
fig_post_KDE = figure('Position',[100,100,400,200],'Color','w');
fig_post_KDE = sampler.PlotPrior(fig_post_KDE,1);
fig_post_KDE = sampler.PlotKDEstimate(fig_post_KDE,1);

fig_post_KDE = figure('Position',[100,350,400,200],'Color','w');
fig_post_KDE = sampler.PlotPrior(fig_post_KDE,2);
fig_post_KDE = sampler.PlotKDEstimate(fig_post_KDE,2);



%% Testing the Gaussian Model and Prior Mapping for Joint Priors
%{
    Finally, we can test the ability to utilize joint priors and
    appropriately map them to the parameters by utilizing a Gaussian
    generating model and a joint Normal-Inverse-Gamma prior. This is also a
    conjugate prior for a Gaussian model, so we can again compare to an
    analytical solution.
%}

model = GaussianModel(1);   % generating model
mu_true = 120; sig_true = 30;   % true parameters we are estimating

% generating the calibration data
N_data = 20;
model.setParams([mu_true,sig_true^2]);
X_data = model.rnd(N_data);

% setting up the normal-inverse-gamma prior
prior = N_IG_Model();
mu0_prior = 100;    % hyper parameters
v0_prior = 0.1;
mu_sig = 50^2;
sig_sig = 30^2;
alpha_prior = (mu_sig/sig_sig)^2 + 1;
beta_prior = mu_sig * (alpha_prior - 1);
prior.setParams([mu0_prior,v0_prior,alpha_prior,beta_prior])
X_prior = prior.random(5000);   % samples from the prior

% setting up the conjugate posterior
mu0_post = (v0_prior*mu0_prior + N_data*mean(X_data))/(v0_prior + N_data);
v0_post = v0_prior + N_data;
alpha_post = alpha_prior + N_data/2;
beta_post = beta_prior + 0.5*sum((X_data-mean(X_data)).^2) + 0.5*N_data*v0_prior*(mean(X_data) - mu0_prior)^2 / (v0_prior + N_data);
posterior = N_IG_Model();
posterior.setParams([mu0_post,v0_post,alpha_post,beta_post]);
X_post = posterior.random(5000);


gammas = linspace(0,1.,10)';   % gamma values to use for the SMC algorithm
transforms = { {2,LogTransform()} };

mapping = {1,[1,2];};

sampler = SMC_Sampler(2000,model,{prior},2.38,X_data,gammas,{'\mu','\sigma'},mapping,transforms);
sampler.GeneratePosterior();

% plotting the convergence of the parameters across the intermediate
% distributions
figs = sampler.PlotConvergence();
figure(figs{1}); hold all
yline(mu_true,'--k','LineWidth',2)
yline(mean(X_post(:,1)),'--b','LineWidth',2)
yline(quantile(X_post(:,1),0.025),':b','LineWidth',2)
yline(quantile(X_post(:,1),0.975),':b','LineWidth',2)

figure(figs{2}); hold all
yline(sig_true.^2,'--k','LineWidth',2)
yline(mean(X_post(:,2)),'--b','LineWidth',2)
yline(quantile(X_post(:,2),0.025),':b','LineWidth',2)
yline(quantile(X_post(:,2),0.975),':b','LineWidth',2)

% Comparing the simulated posterior to the conjugate prior
figs = sampler.PlotPosteriorHistograms();
for i=1:2
figure(figs{i}); hold all
histogram(X_prior(:,i),'FaceColor','r','Normalization','pdf')
histogram(X_post(:,i),'FaceColor','g','Normalization','pdf')
end
xline(sig_true^2,'--k')

