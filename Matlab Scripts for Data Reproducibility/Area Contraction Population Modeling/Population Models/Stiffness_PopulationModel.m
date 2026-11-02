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

        % log likelihood that tries to fit values to the mean
        function logL = log_likelihood(obj, X)
            %{
                This function approximates the likelihood function on the
                summary metrics of the data matrix X. Here X is a Nx3
                matrix, where the first column is the temperature group (-1
                for t=0 group, belonging to all temperature groups), the
                second column is the cross-linking time, and the third is
                the stiffness value. The summary metrics at each time point
            %}

            % parameterizing the sampler and initializing the log
            % likelihood
            M_repet = 200;
            logL = 0;
            use_custom = 1;

            if use_custom
                log_L_norm = @(data, mu, cov) -0.5*( length(mu)*log(2*pi) + log(det(cov)) + (data-mu)' * (cov \ (data-mu)) );
            else
                log_L_norm = @(data, mu, cov) log( mvnpdf(data,mu,cov) );
            end

            %%%%%%% ~~~~~ Step 1: initial data -- direct comparison to the initial data ~~~~~ %%%%%%%
            data_0 = X(X(:,1)==-1,3);        
            mu_E0 = mean(data_0); sig_E0 = std(data_0);
            metrics_E0 = [mu_E0; sig_E0];

            sim_data_0 = obj.rnd_trajectories(0,length(data_0)*M_repet,0);
            sim_data_0 = reshape(sim_data_0,[M_repet,length(data_0)]);
            mu_E0_sim = mean(sim_data_0,2)'; sig_E0_sim = std(sim_data_0,[],2)';
            metrics_E0_sim = [mu_E0_sim; sig_E0_sim];

            mu_metric_E0_sim = mean(metrics_E0_sim,2);
            cov_metric_E0_sim = diag(var(metrics_E0_sim,[],2)); %cov(metrics_E0_sim');
            logL = logL + log_L_norm(metrics_E0, mu_metric_E0_sim, cov_metric_E0_sim);


            %%%%%%% ~~~~~ Step 2: mean time data -- compare fit parameters ~~~~~ %%%%%%%
            X_ = X(X(:,1)>=0,:);    % get all of the except the t=0 data
            temp_groups = unique(X_(:,1));  % get the unique temperature groups

            % function of the mean trajectory, penalizing overly small taus (large 1/tau)
            ft = @(x,t_) x(1).*(1 - exp(-t_./x(2))) + 0*(1e-4)*x(2); 
            % opts = optimset('Display','off');
            % param_fit = @(t_,y_) lsqcurvefit(ft, [60,0.1],t_,y_,[0,0],[],[],[],[],[],[],opts);

            for i=1:length(temp_groups)
                ind_i = find(X_(:,1)==temp_groups(i));
                t_i = X_(ind_i,2); % converting to milliminutes to avoid numerical round off

                % sorting things so that all of the same time groups are together
                [t_i,sorted_inds] = sort(t_i);
                time_groups = unique(t_i);                  % find the unique time values
                N_samples_per_group = sum(t_i==time_groups(i));

                data_i = X_(ind_i,3) - mu_E0;
                data_i = data_i(sorted_inds);
                mean_data = zeros(length(time_groups)+1,1);
                for k=1:length(time_groups)
                    mean_data(k+1) = median(data_i(t_i==time_groups(k)));
                end
                [mu_dY, mu_tau] = obj.FitRelaxation([0;time_groups], mean_data);

                % Fit model to data.
                % fitresult_true = param_fit(t_i,data_i);
                % mu_dY = fitresult(1); mu_inv_tau = fitresult(2); % pull out the model parameters
                sqe_Y = (data_i - ft([mu_dY,mu_tau], t_i)).^2;          % pull out the squared error between data and mean model fit

                data_sim = obj.rnd_trajectories(time_groups',M_repet*N_samples_per_group,temp_groups(i));   % initializing matrix for simulation data
                data_mean_sim = zeros(M_repet,length(time_groups)+1);
                for k=1:M_repet
                    data_mean_sim(k,2:end) = median(data_sim(1+N_samples_per_group*(k-1):N_samples_per_group*k,:),1) - mu_E0_sim(k);
                end

                [mu_dY_sim, mu_tau_sim] = obj.FitRelaxation([0;time_groups],data_mean_sim');

                std_Y_sim = zeros(M_repet,length(time_groups));
                for k=1:M_repet
                    std_Y_sim(k,:) = sqrt(mean( (data_sim(1+N_samples_per_group*(k-1):N_samples_per_group*k,:) - mu_E0_sim(k) - ft([mu_dY_sim(k),mu_tau_sim(k)],time_groups)').^2, 1) );
                end

                std_Y = zeros(length(time_groups),1);       % preallocating space for data variance at the different groups
                j0 = 1;
                for j=1:length(time_groups)
                    j1 = j0+sum(t_i==time_groups(j))-1;
                    std_Y(j) = sqrt(mean(sqe_Y(j0:j1)));  % calculate the variance of the data (compared to the model fit) at each time point
                    j0=j1+1;
                end

                % setting up the metrics
                curve_metrics = [mu_dY; mu_tau; std_Y];
                sim_metrics = [mu_dY_sim; mu_tau_sim; std_Y_sim'];

                mu_metric_curve_sim = mean(sim_metrics,2);  % mean of simulated metrics
                % cov_metric_curve_sim = cov(sim_metrics');   % emperical covariance of the simulated metrics
                % cov_metric_curve_sim = cov_metric_curve_sim + (1e-8)*diag(ones(size(mu_metric_curve_sim))); % regularizing
                cov_metric_curve_sim = diag(var(sim_metrics,[],2) + (1e-5)*sum(var(sim_metrics,[],2)));
                logL = logL + log_L_norm(curve_metrics, mu_metric_curve_sim, cov_metric_curve_sim);

            end
        end

        function [dY, tau] = FitRelaxation(obj,t,Y)
            %{
                Uses an integral form of the solution y(t) =
                dY*(1-exp(-t/tau)) to estimate tau and dY. Y can be a
                matrix where each row corresponds to a single sample of
                y(t)
            %}
            m = length(t);  % number of time points
            k = size(Y,2);  % number of samples being fit

            I_Y = cumtrapz(t,Y,1);
            A = zeros(m,2,k);
            A(:,1,:) = I_Y;
            A(:,2,:) = repmat(t,[1,1,k]);
            A_sq = pagemtimes(pagetranspose(A),A);
            Y = reshape(Y,[size(Y,1),1,size(Y,2)]);
            p_ = pagemldivide(A_sq,pagemtimes(pagetranspose(A), Y));

            tau = -1./p_(1,:);
            dY = -p_(2,:)./p_(1,:);
        end

        % function logL = log_likelihood(obj, X)
        %     %{
        %         This function approximates the likelihood function on the
        %         summary metrics of the data matrix X. Here X is a Nx3
        %         matrix, where the first column is the temperature group (-1
        %         for t=0 group, belonging to all temperature groups), the
        %         second column is the cross-linking time, and the third is
        %         the stiffness value. The summary metrics at each time point
        %         is the group mean and standard deviation.
        %     %}
        % 
        %     M_repet = 100; % use M_repet datasets sampled from the model to approximate the likelihood
        %     logL = 0; % initializing
        % 
        %     % breaking out the temperature groups and the time groups
        %     temp_groups = unique(X(:,1));
        %     for i=1:length(temp_groups)
        %         % looping through all of the temperature groups
        %         T_i = max(temp_groups(i),0); % sets the t=0 group to be from the control group
        %         X_i = X(X(:,1)==temp_groups(i),:); % getting just the data for this temperature group
        %         time_groups = unique(X_i(:,2));
        % 
        % 
        %         for j=1:length(time_groups)
        %             % looping through all of the time groups
        %             t_j = time_groups(j);
        %             X_j = X_i(X_i(:,2)==t_j,3); % pull out just the stiffness values
        %             metrics_data = [median(X_j); iqr(X_j)];
        % 
        %             X_j_sim = obj.rnd_trajectories(t_j,M_repet*length(X_j),T_i); % generate all of the synthetic data
        %             X_j_sim = reshape(X_j_sim,[length(X_j),M_repet]);
        %             metrics_sim = [median(X_j_sim,1); iqr(X_j_sim,1)]; % calculate the metrics for all repetitions
        %             mu_metrics = mean(metrics_sim,2); % calculate the mean and covariance for the simulation metrics
        %             cov_metrics = cov(metrics_sim');
        % 
        %             logL = logL + log(mvnpdf(metrics_data, mu_metrics, cov_metrics)); % update the total log likelihood for the new time point
        %         end
        %     end
        % end


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

            % X_data is the 3 column matrix of [temp group, time,
            % stiffness]. temp group == -1 is the initial stiffness data

            passes = zeros(M_repet,1);

            % Step 1: initial data -- direct comparison to the initial data
            ind_0 = find(X(:,1)==-1);
            data_0 = X(ind_0,3);        
            mu_E0 = mean(data_0);
            std_E0 = std(data_0);

            sim_data_0 = obj.rnd_trajectories(0,length(data_0)*M_repet,0);
            sim_data_0 = reshape(sim_data_0,[M_repet,length(data_0)]);
            mu_E0_sim = mean(sim_data_0,2);
            std_E0_sim = std(sim_data_0,[],2);
            
            % Are these simulations close enough?
            passes = passes + ( abs((mu_E0 - mu_E0_sim)./mu_E0) < errors(1));
            passes = passes + ( abs((std_E0 - std_E0_sim)./std_E0) < errors(2));

            % Step 2: mean time data -- compare fit parameters
            X_ = X(X(:,1)>=0,:);
            temp_groups = unique(X_(:,1));
            
            % ft = fittype( 'a*(1-exp(-b*x))', 'independent', 'x', 'dependent', 'y' );
            % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            % opts.Display = 'Off';
            % opts.Lower = [0 0];
            % opts.StartPoint = [60 0.1];
            ft = @(x,t_) x(1).*(1 - exp(-x(2).*t_));

            opts = optimset('Display','off');
            param_fit = @(t_,y_) lsqcurvefit(ft, [60,0.1],t_,y_,[0,0],[],[],[],[],[],[],opts);

            for i=1:length(temp_groups)
                ind_i = find(X_(:,1)==temp_groups(i));
                t_i = X_(ind_i,2);
                
                % sorting things so that all of the same time groups are
                % together
                [t_i,sorted_inds] = sort(t_i);

                data_i = X_(ind_i,3) - mu_E0;
                data_i = data_i(sorted_inds);
               
                % Fit model to data.
                % fitresult = fit( t_i, data_i, ft, opts );
                fitresult = param_fit(t_i,data_i);
                mu_dY = fitresult(1); mu_inv_tau = fitresult(2); % pull out the model parameters
                sqe_Y = (data_i - ft(fitresult,t_i)).^2;          % pull out the squared error between data and mean model fit

                time_groups = unique(t_i);                  % find the unique time values
                data_sim = zeros(M_repet,length(data_i));   % initializing matrix for simulation data
                var_Y = zeros(length(time_groups),1);       % preallocating space for data variance at the different groups

                j0 = 1;
                for j=1:length(time_groups)
                    j1 = j0+sum(t_i==time_groups(j))-1;
                    sim_data_i = obj.rnd_trajectories(time_groups(j),sum(t_i==time_groups(j))*M_repet,temp_groups(i));       % generate the data
                    data_sim(:,j0:j1) = reshape(sim_data_i,[M_repet,sum(t_i==time_groups(j))]) - mu_E0_sim;     % reshape to match the experimental data
                    
                    var_Y(j) = mean(sqe_Y(j0:j1));  % calculate the variance of the data (compared to the model fit) at each time point

                    j0=j1+1;
                end

                mu_dY_sim = zeros(M_repet,1); mu_inv_tau_sim = zeros(M_repet,1);
                var_Y_sim = zeros(M_repet,length(time_groups));
                for k=1:M_repet
                    try
                        % fitresult = fit( t_i, data_sim(k,:)', ft, opts );
                        fitresult = param_fit( t_i, data_sim(k,:)');
                    catch
                        disp("Something fucky is happening...")
                    end
                    mu_dY_sim(k) = fitresult(1); mu_inv_tau_sim(k) = fitresult(2);
                    sqe_Y_sim = ( data_sim(k,:) - ft(fitresult,t_i)' ).^2;
                    j0 = 1;
                    for j=1:length(time_groups)
                        j1 = j0+sum(t_i==time_groups(j))-1;
                        var_Y_sim(k,j) = mean(sqe_Y_sim(j0:j1));
                        j0=j1;
                    end
                end

                % Are these simulations close enough?
                passes = passes + ( abs((mu_dY - mu_dY_sim)./mu_dY) < errors(3));
                passes = passes + ( abs((mu_inv_tau - mu_inv_tau_sim)./mu_inv_tau) < errors(4));
                passes = passes + ( mean( abs((var_Y' - var_Y_sim)./var_Y'), 2) < errors(5));
            end
            sims_pass = passes == (2 + 3*length(temp_groups)); % did at least half of the metrics pass?

            approx_L = sum(sims_pass) / M_repet;

            % 
            % 
            % for i=1:N_experiments
            %     % for each of the time points
            %     N_samps = X{i}{1};      % simulation sample size
            %     Sim_params = X{i}{2};   % cross_linking time
            %     X_data = X{i}{3};       % stiffness data
            % 
            %     E_data = X_data(:,1);              % Stiffness
            %     metrics_data = [mean(E_data,1), std(E_data,1)];
            % 
            %     % generate all of the trajectories for the replicates
            %     N_total = M_repet * N_samps;
            %     sim_data = obj.rnd_trajectories([Sim_params(1)],N_total,Sim_params(2));
            %     E_sim = sim_data(:,1);             % Stiffness
            % 
            %                     % metrics of interest: median and IQR of initial and final
            %     % mass
            %     means = zeros(M_repet,1);   
            %     stds = zeros(M_repet,1);
            % 
            %     for j=1:M_repet
            %         means(j,:) = mean(E_sim(1+N_samps*(j-1):N_samps*j,:),1);
            %         stds(j,:)    = std(E_sim(1+N_samps*(j-1):N_samps*j,:),1);
            %     end
            %     metrics = [means'; stds']'; % median M0, median dM, IQR M0, IQR dM
            %     err_ = (metrics - metrics_data);
            % 
            %     % normalize both M0 and dM metrics by the M0 data to prevent divide by 0 if dM=0
            %     err_(:,1) = err_(:,1) ./ metrics_data(:,1); 
            %     err_(:,2) = err_(:,2) ./ metrics_data(:,2); 
            % 
            %     err_central(:,i) = err_(:,1);
            %     err_spread(:,i) = err_(:,2);
            % 
            %     % experiments_pass = experiments_pass .* prod((abs(err_) < errors),2);
            % end
            % 
            % % RMSE error across all four time points
            % total_err_central = 100*sqrt(mean(err_central.^2,2));
            % total_err_spread = 100*sqrt(mean(err_spread.^2,2));
            % experiments_pass = (total_err_central < errors(1)).*(total_err_spread < errors(2));
            % 
            % approx_L = sum(experiments_pass) / M_repet;

        end

    end


end