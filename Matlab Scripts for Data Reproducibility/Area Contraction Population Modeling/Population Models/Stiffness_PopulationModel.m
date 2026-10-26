%% Population Model for Time-Varying Stiffness of the Hydrogel
%{
    Last Updated: MJB 10/17/2026

    This specific class inherits from the PopulationModel class and
    implements the approximate likelihood for the stiffness simulations 
    needed to estimate the stiffness parameters using SMC ABC
%}

classdef Stiffness_PopulationModel < PopulationModel
    %{

    %}

    properties

    end

    methods
        function obj = Stiffness_PopulationModel(include_effect)
            obj@PopulationModel(include_effect)
        end

        function approx_L = approx_likelihood(obj,X,errors,M_repet)
            %{
                This function returns the approximate likelihood of the
                mass population model given the initial mass and mass at a
                given cross-linking time
                
                The dataset (X) includes data from the 0, 2, 4, and 6 hour
                time groups as individual "experiments." Therefore, X has 4
                fours, where each row has {sample_size, [cross_linking_time, 
                temperature_group], [initial_mass, final_mass]}
            %}

            N_experiments = size(X,1);
            err_central = zeros(M_repet,N_experiments);
            err_spread = zeros(M_repet,N_experiments);

            for i=1:N_experiments
                % for each of the time points
                N_samps = X{i}{1};      % simulation sample size
                Sim_params = X{i}{2};   % cross_linking time
                X_data = X{i}{3};       % stiffness data

                E_data = X_data(:,1);              % Stiffness
                metrics_data = [mean(E_data,1), std(E_data,1)];

                % generate all of the trajectories for the replicates
                N_total = M_repet * N_samps;
                sim_data = obj.rnd_trajectories([Sim_params(1)],N_total,Sim_params(2));
                E_sim = sim_data(:,1);             % Stiffness

                                % metrics of interest: median and IQR of initial and final
                % mass
                means = zeros(M_repet,1);   
                stds = zeros(M_repet,1);

                for j=1:M_repet
                    means(j,:) = mean(E_sim(1+N_samps*(j-1):N_samps*j,:),1);
                    stds(j,:)    = std(E_sim(1+N_samps*(j-1):N_samps*j,:),1);
                end
                metrics = [means'; stds']'; % median M0, median dM, IQR M0, IQR dM
                err_ = (metrics - metrics_data);

                % normalize both M0 and dM metrics by the M0 data to prevent divide by 0 if dM=0
                err_(:,1) = err_(:,1) ./ metrics_data(:,1); 
                err_(:,2) = err_(:,2) ./ metrics_data(:,2); 

                err_central(:,i) = err_(:,1);
                err_spread(:,i) = err_(:,2);

                % experiments_pass = experiments_pass .* prod((abs(err_) < errors),2);
            end
            
            % RMSE error across all four time points
            total_err_central = 100*sqrt(mean(err_central.^2,2));
            total_err_spread = 100*sqrt(mean(err_spread.^2,2));
            experiments_pass = (total_err_central < errors(1)).*(total_err_spread < errors(2));

            approx_L = sum(experiments_pass) / M_repet;

        end

    end


end