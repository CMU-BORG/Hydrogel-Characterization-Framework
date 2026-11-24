%% Population Model for Time-Varying Stiffness of the Hydrogel
%{
    Last Updated: MJB 11/7/2026

    This specific class inherits from the PopulationModel class and
    implements the approximate likelihood for the stiffness simulations 
    needed to estimate the stiffness parameters using SMC ABC
%}

classdef Stiffness_PopulationModel < PopulationModel
    %{
        Specific implementation of the population model for mass
        measurements, with the likelihood implemented with as a series of
        kernel density estimates at each time and temperature setting.
    %}

    properties

    end

    methods
        function obj = Stiffness_PopulationModel(include_effect)
            obj@PopulationModel(include_effect)
        end

        function logL = log_likelihood(obj,X)
            %{
                Synethic likelihood estimate using kernel density estimates
                at each time and temperature group.
            %}
            M_samples_per_KDE = 1000;   % number of samples generated to construct the kernel density estimate 

            X(X(:,1)<0,1) = 0;  % re-assigning the control group to the 27C group. Does not change anything because t=0 response is temperature independent

            temp_groups = unique(X(:,1)); % find all of the temperature groups to check

            logL = 0;   % initializing the log-likelihood
            log_norm = @(x,samps) log( ksdensity(samps,x) );    % log kernel density pdf estimate of the points x given the samples samps

            for j=1:length(temp_groups)
                % loop through all of the temperature groups

                data_j = X(X(:,1)==temp_groups(j),:);   % pull out the data associated with the current temperature group
                time_groups = unique(data_j(:,2));      % find all of the times in the current temperature group

                for m=1:length(time_groups)
                    % loop through all of the time groups

                    data_m = data_j(data_j(:,2)==time_groups(m),3:end);     % pull out the data associated with the current time group

                    % generate the simulated data for kernel density
                    % estimate. For stiffness data, we need only the t=t_i 
                    % points for each sample
                    sim_data = obj.rnd_trajectories(time_groups(m),M_samples_per_KDE,temp_groups(j)); 
                    
                    logL = logL + sum(log_norm(data_m, sim_data));      % calculate the log kernel density estimate and update the log likelihood
                end
            end
        end

    end


end