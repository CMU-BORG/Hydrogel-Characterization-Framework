%% Testing SMC Sampler 
%{

%}

clear all; close all; clc;

% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool("local", 8);  % Adjust the number of workers if needed
end

timestamp = datestr(now, 'yyyymmdd_HHMM');

%% Testing Using Bernoulli Model and Conjugate Prior
bernoulli_model = BernoulliModel(0);


N_data = 20;    % number of datapoints to use in the dataset
p_true = 0.75;  % true value of the model parameter for generating the data
bernoulli_model.setParams([p_true]);
X_data_true = bernoulli_model.rnd(N_data);

figure; histogram(X_data_true,'Normalization','probability','BinWidth',0.333);

a_prior = 1; b_prior = 1; % hyper parameters of the prior
priors = {makedist("Beta","a",a_prior,"b",b_prior)}; % using a Beta prior so we can test the SMC algorithm against a known conjugate prior
posterior = makedist("Beta","a",a_prior+sum(X_data_true),"b",b_prior+N_data-sum(X_data_true));    % true posterior based on conjugacy

gammas = min(linspace(0,1.,10),1)';   % gamma values to use for the SMC algorithm
transforms = {{1, LogitTransform(0,1)};};

sampler = SMC_Sampler(1000,bernoulli_model,priors,2.38,X_data_true,gammas,["p"],{},transforms);
sampler.GeneratePosterior();

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
drawnow()

%%
% Comparing to the kernel density estimates of the PDFs
fig_post_KDE = figure('Position',[100,100,400,200],'Color','w');
fig_post_KDE = sampler.PlotPrior(fig_post_KDE,1);
fig_post_KDE = sampler.PlotKDEstimate(fig_post_KDE,1);




%% Testing the ABC SMC method using the Bernoulli Model

data_1 = X_data_true(1:N_data/2);
data_2 = X_data_true(1+N_data/2:end);

X_data = {{N_data/2,[],sum(data_1)};
          {N_data/2,[],N_data/2 - sum(data_2)}};    % data cell array for ABC method

err_ = linspace(0.5^(1/4),0.005^(1/4),15).^4;
gammas = (N_data/2)*[err_',err_'];   % gamma values to use for the SMC algorithm


sampler = ABC_SMC_Sampler(1000,10,bernoulli_model,priors,2.38,X_data,gammas,["p"],{},transforms);
sampler.GeneratePosterior();

figs = sampler.PlotConvergence();
figure(figs{1}); hold all
yline(p_true,'--k','LineWidth',2)
yline(posterior.mean(),'-b','LineWidth',2)
yline(posterior.icdf([0.025,0.975]),':b','LineWidth',2)

fig_post = sampler.PlotPosteriorHistograms();
x_ = linspace(0,1,300);
figure(fig_post{1}); hold all
plot(x_,priors{1}.pdf(x_),'--','Color','r','LineWidth',2)
plot(x_,posterior.pdf(x_),'-','Color','b','LineWidth',2);
plot(x_,sampler.pdf(x_'),':g','LineWidth',2)
drawnow()

%% Testing with Gamma for multi-variable fitting

gamma_model = GammaModel(0);

N_data = 30;    % number of datapoints to use in the dataset
alpha_true = 1;  % true value of the model parameters for generating the data
beta_true = 10;
gamma_model.setParams([alpha_true,beta_true]);
X_data = gamma_model.rnd(N_data);

figure; histogram(X_data,'Normalization','probability');

alpha_prior_a = 0; alpha_prior_b = 5; % hyper parameters of the prior
beta_prior_a = 0; beta_prior_b = 50; 

% numerically calculating the true posterior
N_points = 100;
alpha_test = linspace(alpha_prior_a,alpha_prior_b,N_points+1);
alpha_test(1) = []; % removing 0
beta_test = linspace(beta_prior_a,beta_prior_b,N_points+1);
beta_test(1) = [];

[A,B] = meshgrid(alpha_test,beta_test);
post_ = 0*A;
post_approx = 0*A;
for i=1:N_points
    for j=1:N_points
        gamma_model.setParams([A(i,j),B(i,j)]);
        post_(i,j) = exp(gamma_model.log_likelihood(X_data)); % no prior mult because it is just constant
        post_approx(i,j) = gamma_model.approx_likelihood({{N_data, [], [mean(X_data),std(X_data)]};},[0.05,0.05],20);
    end
end

I = trapz(beta_test, trapz(alpha_test, post_, 2));
I_approx = trapz(beta_test, trapz(alpha_test, post_approx, 2));


post_ = post_ / I;
post_approx = post_approx / I_approx;

f_alpha = trapz(beta_test, post_, 1);
f_beta = trapz(alpha_test, post_, 2);
F_alpha = cumtrapz(alpha_test,f_alpha); % CDF
F_beta = cumtrapz(beta_test,f_beta); % CDF

val_Fa = (F_alpha>0.001)&(F_alpha<0.999);
val_Fb = (F_beta>0)&(F_beta<1);
alpha_inv_CDF = @(p) interp1(F_alpha(val_Fa),alpha_test(val_Fa),p);
beta_inv_CDF = @(p) interp1(F_beta(val_Fb),beta_test(val_Fb),p);
alpha_CI = alpha_inv_CDF([0.025,0.975]);
beta_CI = beta_inv_CDF([0.025,0.975]);
alpha_med = alpha_inv_CDF(0.5);
beta_med = beta_inv_CDF(0.5);


f_alpha_approx = trapz(beta_test, post_approx, 1);
f_beta_approx = trapz(alpha_test, post_approx, 2);

figure; hold all
plot(alpha_test,f_alpha,'-r','LineWidth',2);
plot(alpha_test,f_alpha_approx,'--r','LineWidth',2);
plot(beta_test,f_beta,'-b','LineWidth',2);
plot(beta_test,f_beta_approx,'--b','LineWidth',2);


priors = {makedist("Uniform","lower",alpha_prior_a,"upper",alpha_prior_b),
          makedist("Uniform","lower",beta_prior_a,"upper",beta_prior_b)}; % using a Beta prior so we can test the SMC algorithm against a known conjugate prior

proposals = [makedist("Normal","mu",0,"sigma",0.5),
             makedist("Normal","mu",0,"sigma",5)];
gammas = min(linspace(0,1.,10),1)';   % gamma values to use for the SMC algorithm
transforms = cell(2,1);
transforms{1} = {1,LogTransform()};
transforms{2} = {2,LogTransform()};

sampler = SMC_Sampler(1000,gamma_model,priors,2.38,X_data,gammas,["\alpha","\beta"],{},transforms);
sampler.GeneratePosterior();

%%

% f_post = sampler.pdf([A(:),B(:)]);
% f_post = reshape(f_post,size(A));
% fa_post = trapz(beta_test,f_post,1);
% fb_post = trapz(alpha_test,f_post,2);
fa_post = sampler.marginal_pdf(1,alpha_test');
fb_post = sampler.marginal_pdf(2,beta_test');

mu_data = mean(X_data);
sig_data = std(X_data);
alpha_data = (mu_data/sig_data)^2;
beta_data = sig_data^2 / mu_data;

% figs = sampler.PlotConvergence();
% figure(figs{1}); hold all
% yline(alpha_true,'--k','LineWidth',2)
% yline(alpha_data,'--g','LineWidth',2);
% yline(alpha_med,'--b','LineWidth',2);
% yline(alpha_CI,':b','LineWidth',2)
% 
% figure(figs{2}); hold all
% yline(beta_true,'--k','LineWidth',2)
% yline(beta_data,'--g','LineWidth',2);
% yline(beta_med,'--b','LineWidth',2);
% yline(beta_CI,':b','LineWidth',2)

post_resamps = sampler.random(5000);

fig_post = sampler.PlotPosteriorHistograms();
x_ = linspace(0,50,300);
figure(fig_post{1}); hold all
xline(alpha_true,'--k','LineWidth',2);
xline(alpha_data,'--g','LineWidth',2);
plot(0.05*x_,priors{1}.pdf(0.05*x_),'--','Color','r','LineWidth',2)
plot(alpha_test,f_alpha,'--b','LineWidth',2)
plot(alpha_test,fa_post,':g','LineWidth',2)
histogram(post_resamps(:,1),'FaceAlpha',0.2,'Normalization','pdf')

figure(fig_post{2}); hold all
xline(beta_true,'--k','LineWidth',2);
xline(beta_data,'--g','LineWidth',2);
plot(x_,priors{2}.pdf(x_),'--','Color','r','LineWidth',2)
plot(beta_test,f_beta,'--b','LineWidth',2)
plot(beta_test,fb_post,':g','LineWidth',2)
histogram(post_resamps(:,2),'FaceAlpha',0.2,'Normalization','pdf')


% sampler.PlotJointPosterior();
% nexttile(3); hold all
% contour(A,B,post_,15)

%% % Comparing to the kernel density estimates of the PDFs
fig_post_KDE = figure('Position',[100,100,400,200],'Color','w');
fig_post_KDE = sampler.PlotPrior(fig_post_KDE,1);
fig_post_KDE = sampler.PlotKDEstimate(fig_post_KDE,1);

fig_post_KDE = figure('Position',[100,350,400,200],'Color','w');
fig_post_KDE = sampler.PlotPrior(fig_post_KDE,2);
fig_post_KDE = sampler.PlotKDEstimate(fig_post_KDE,2);


%% Testing ABC SMC method using the Gamma Model
% using the same model data from previous gamma model so we can use the
% same numerically calculated posterior to comparison

gamma_model_approx = GammaModel(0,1);
sampler = SMC_Sampler(1000,gamma_model_approx,priors,2.38,X_data,gammas,["\alpha","\beta"],{},transforms);
sampler.GeneratePosterior();


figs = sampler.PlotConvergence();
figure(figs{1}); hold all
yline(alpha_true,'--k','LineWidth',2)
yline(alpha_data,'--g','LineWidth',2);
yline(alpha_med,'--b','LineWidth',2);
yline(alpha_CI,':b','LineWidth',2)

figure(figs{2}); hold all
yline(beta_true,'--k','LineWidth',2)
yline(beta_data,'--g','LineWidth',2);
yline(beta_med,'--b','LineWidth',2);
yline(beta_CI,':b','LineWidth',2)

fig_post = sampler.PlotPosteriorHistograms();
x_ = linspace(0,50,300);
figure(fig_post{1}); hold all

xline(alpha_true,'--k','LineWidth',2);
xline(alpha_data,'--g','LineWidth',2);
plot(0.05*x_,priors{1}.pdf(0.05*x_),'--','Color','r','LineWidth',2)
plot(alpha_test,f_alpha,'--b','LineWidth',2)


figure(fig_post{2}); hold all
xline(beta_true,'--k','LineWidth',2);
xline(beta_data,'--g','LineWidth',2);
plot(x_,priors{2}.pdf(x_),'--','Color','r','LineWidth',2)
plot(beta_test,f_beta,'--b','LineWidth',2)

for i=1:2
    figure(fig_post{i}); hold all
    histogram(priors{i}.random(1000,1),'Normalization','pdf')
end


%% Testing the Gaussian Model and Prior Mapping for Joint Priors

model = GaussianModel(1);
mu_true = 120;
sig_true = 30;

N_data = 6;
model.setParams([mu_true,sig_true^2]);
X_data = model.rnd(N_data);

prior = N_IG_Model();
mu0_prior = 100;
v0_prior = 0.1;
mu_sig = 50^2;
sig_sig = 30^2;
alpha_prior = (mu_sig/sig_sig)^2 + 1;
beta_prior = mu_sig * (alpha_prior - 1);
prior.setParams([mu0_prior,v0_prior,alpha_prior,beta_prior])
mapping = {1,[1,2];};
X_prior = prior.random(5000);

figure;
subplot(1,2,1); histogram(X_prior(:,1),"Normalization","pdf");
subplot(1,2,2); histogram(X_prior(:,2),"Normalization","pdf");

% conjugate posterior
mu0_post = (v0_prior*mu0_prior + N_data*mean(X_data))/(v0_prior + N_data);
v0_post = v0_prior + N_data;
alpha_post = alpha_prior + N_data/2;
beta_post = beta_prior + 0.5*sum((X_data-mean(X_data)).^2) + 0.5*N_data*v0_prior*(mean(X_data) - mu0_prior)^2 / (v0_prior + N_data);
posterior = N_IG_Model();
posterior.setParams([mu0_post,v0_post,alpha_post,beta_post]);

gammas = min(linspace(0,1.,10),1)';   % gamma values to use for the SMC algorithm
transforms = cell(1,1);
transforms{1} = {2,LogTransform()};

sampler = SMC_Sampler(2000,model,{prior},2.38,X_data,gammas,{'\mu','\sigma'},mapping,transforms);
sampler.GeneratePosterior();

X_post = posterior.random(5000);

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

figs = sampler.PlotPosteriorHistograms();
for i=1:2
figure(figs{i}); hold all
histogram(X_prior(:,i),'FaceColor','r','Normalization','pdf')
histogram(X_post(:,i),'FaceColor','g','Normalization','pdf')
end
xline(sig_true^2,'--k')

%% Gaussian with ABC SMC

metrics = [mean(X_data), std(X_data,1)];
X_data_cell = {{N_data, [], metrics};};

N_gens = 10;
n_ = 3;
gammas_mu = linspace(1,0.025^(1/n_),N_gens).^n_;
gammas_sig = gammas_mu; % linspace(4,0.01,N_gens);
gammas = zeros(N_gens,1,2);
gammas(:,1,1) = gammas_mu;
gammas(:,1,2) = gammas_sig;

sampler = ABC_SMC_Sampler(2000,5,model,{prior},2.38,X_data_cell,gammas,["\mu","\sigma"],mapping,transforms);
sampler.GeneratePosterior();

figs = sampler.PlotConvergence();
figure(figs{1}); hold all
yline(mu_true,'--k','LineWidth',2)
yline(mu_true,'--k','LineWidth',2)
yline(mean(X_post(:,1)),'--b','LineWidth',2)
yline(quantile(X_post(:,1),0.025),':b','LineWidth',2)
yline(quantile(X_post(:,1),0.975),':b','LineWidth',2)

figure(figs{2}); hold all
yline(sig_true.^2,'--k','LineWidth',2)
yline(mean(X_post(:,2)),'--b','LineWidth',2)
yline(quantile(X_post(:,2),0.025),':b','LineWidth',2)
yline(quantile(X_post(:,2),0.975),':b','LineWidth',2)

fig_post = sampler.PlotPosteriorHistograms();
for i=1:2
figure(fig_post{i}); hold all
histogram(X_prior(:,i),'FaceColor','r','Normalization','pdf')
histogram(X_post(:,i),'FaceColor','g','Normalization','pdf')
end
xline(sig_true^2,'--k')
