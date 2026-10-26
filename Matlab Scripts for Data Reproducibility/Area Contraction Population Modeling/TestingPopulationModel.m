%% Testing Population Model
%{

%}

clear all; close all; clc

pop_model = PopulationModel();
pop_model.set_mu_P(200,50);
pop_model.set_cov_P(20,8,0);
pop_model.set_tau_params(120,30);
pop_model.set_P_effect(0.20,0.05);
pop_model.set_tau_effect(0.20,0.05);

theta_1 = pop_model.rnd_params(1000,0);
theta_2 = pop_model.rnd_params(1000,10);
params = ["P_0","\Delta P","\tau"];

figure("Position",[100,100,1200,400],'Color','w');
for i=1:3
    subplot(1,3,i); hold all
    histogram(theta_1(:,i),25,'Normalization','pdf');
    histogram(theta_2(:,i),25,'Normalization','pdf');
    xlabel(params(i))
    set(gca,"FontSize",14,"FontName","Calibri")
end

t = linspace(0,1000,100);
P_1 = pop_model.rnd_trajectories(t,100,0);
P_2 = pop_model.rnd_trajectories(t,100,10);

figure("Position",[100,100,1200,400],'Color','w');
for i=1:100
    subplot(1,2,1); hold all
    plot(t,P_1(i,:),'-','Color',[1,0.25,0.25]*rand(1),'LineWidth',1);
    
    subplot(1,2,2); hold all
    plot(t,P_2(i,:),'-','Color',[0.25,0.25,1]*rand(1),'LineWidth',1);
end
for i=1:2
    subplot(1,2,i);
    ylim([150,400])
    set(gca,"FontSize",14,"FontName","Calibri")
    xlabel("Time [min]")
    ylabel("P(t)")
end


