%% Population Model for Time-Varying Mass of the Hydrogel
%{
    Last Updated: MJB 10/14/2026

    This specific class inherits from the PopulationModel class and
    implements the approximate likelihood for the mass simulations needed
    to estimate the mass parameters using SMC ABC
%}

classdef Mass_PopulationModel < PopulationModel
    %{

    %}

    properties

    end

    methods
        function obj = Mass_PopulationModel(include_effect)
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

            experiments_pass = ones(M_repet,1);
            for i=1:4
                % for each of the time points
                N_samps = X{i}{1};      % simulation sample size
                Sim_params = X{i}{2};   % cross_linking time
                X_data = X{i}{3};       % initial mass and final mass vectors

                M0_data = X_data(:,1);              % initial mass
                dM_data = X_data(:,2) - M0_data;    % change in mass for the given cross-linking time
                M_data = [M0_data, dM_data];
                metrics_data = [median(M_data,1), iqr(M_data,1)];

                % generate all of the trajectories for the replicates
                N_total = M_repet * N_samps;
                sim_data = obj.rnd_trajectories([0,Sim_params(1)],N_total,Sim_params(2));
                M0_sim = sim_data(:,1);             % initial simulated mass
                dM_sim = sim_data(:,2) - M0_sim;    % initial simulated change in mass
                M_sim = [M0_sim, dM_sim];

                % metrics of interest: median and IQR of initial and final
                % mass
                medians = zeros(M_repet,2);   
                IQRs = zeros(M_repet,2);

                for j=1:M_repet
                    medians(j,:) = median(M_sim(1+N_samps*(j-1):N_samps*j,:),1);
                    IQRs(j,:)    = iqr(M_sim(1+N_samps*(j-1):N_samps*j,:),1);
                end
                metrics = [medians'; IQRs']'; % median M0, median dM, IQR M0, IQR dM
                err_ = (metrics - metrics_data);

                % normalize both M0 and dM metrics by the M0 data to prevent divide by 0 if dM=0
                err_(:,1) = err_(:,1) ./ metrics_data(:,1); 
                % err_(:,2) = err_(:,2) ./ metrics_data(:,1); 
                err_(:,3) = err_(:,3) ./ metrics_data(:,3);
                % err_(:,4) = err_(:,4) ./ metrics_data(:,3); 

                avg_error = 0.5*(err_(:,1) + err_(:,2)); %mean(err_,2);
                

                experiments_pass = experiments_pass .* (abs(avg_error) < errors); %prod((abs(err_) < errors),2);
            end
            approx_L = sum(experiments_pass) / M_repet;

        end

    end


end