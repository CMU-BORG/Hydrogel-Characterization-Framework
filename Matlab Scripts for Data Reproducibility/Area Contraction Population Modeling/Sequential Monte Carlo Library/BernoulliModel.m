%% Bernoulli Model class for testing SMC and ABC SMC classes
%{

%}

classdef BernoulliModel < handle
    %{

    %}

    properties
        p       % probability of a 1
        use_approx_likelihood       % flag set to determine if true or approx likelihood should be used
    end

    methods
        function obj = BernoulliModel(use_approx_likelihood)
            obj.p = 0.5; % initializing
            obj.use_approx_likelihood = use_approx_likelihood;
        end

        function setParams(obj,theta)
            % updating the parameter
            obj.p = theta(1);
        end

        function logL = log_likelihood(obj,X)

            if ~obj.use_approx_likelihood
                % log likelihood of the dataset X
                N1 = sum(X);    % number of datapoints equal to 1
                N0 = length(X) - N1;    % number of datapoints equal to 0
                
                logL = N1*log(obj.p) + N0*log(1-obj.p);
            else
                % approximate likelihood 
                
            end
        end

        function X = rnd(obj,N)
            % generates N random Bernoulli variables
            if size(N)==1
                u_ = rand(N,1);
            else
                u_ = rand(N);
            end
            X = u_ < obj.p;
        end

        function approx_L = approx_likelihood(obj,X,errors,M_repet)
            % an approximate likelihood for this model to test the ABC SMC
            % code. In this approximate method, N samples are generated
            % from the model, and the total numbers of 1s are calculated as
            % the summary metric. Then the percent error between the
            % simulation and experimental summary metrics is calculated and
            % compared to the error threshold. This is conducted for 20
            % repetitions

            % M_repet = 500; % number of experiments to run; this can be tuned
                
            experiment_pass = zeros(M_repet,2);

            % Experiment 1: Number of 1s
            N_samps = X{1}{1}; % sample size of the experimental data 
            Sim_params = X{1}{2};  % any input parameters you experiments need; here, we have none
            X_data = X{1}{3};  % summary metrics

            x_ = obj.rnd([N_samps,M_repet]);  % generate data
            N_sim = sum(x_,1);        % number of 1s in the simulated data

            err_i = abs(N_sim - X_data);
            experiment_pass(:,1) = err_i <= errors(1); % is simulated summary metric sufficiently close?

            % Experiment 2: Number of 0s
            N_samps = X{2}{1}; % sample size of the experimental data 
            Sim_params = X{2}{2};  % any input parameters you experiments need; here, we have none
            X_data = X{2}{3};  % summary metrics
            
            x_ = obj.rnd([N_samps,M_repet]);  % generate data
            N_sim = N_samps - sum(x_,1);        % number of 0s in the simulated data

            err_i = abs(N_sim - X_data);
            experiment_pass(:,2) = err_i <= errors(2); % is simulated summary metric sufficiently close?               
            
            all_pass = prod(experiment_pass,2);
            approx_L = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
        end

        function d_approx_L = delta_approx_likelihood(obj,X,errors)
            
            M_repet = 500; % number of experiments to run; this can be tuned
                
            experiment_pass = zeros(M_repet,2);

            % Experiment 1: Number of 1s
            N_samps = X{1}{1}; % sample size of the experimental data 
            Sim_params = X{1}{2};  % any input parameters you experiments need; here, we have none
            X_data = X{1}{3};  % summary metrics

            x_ = obj.rnd([N_samps,M_repet]);  % generate data
            N_1 = sum(x_,1);        % number of 1s in the simulated data
            err_1 = abs(N_1 - X_data);

            % Experiment 2: Number of 0s
            N_samps = X{2}{1}; % sample size of the experimental data 
            Sim_params = X{2}{2};  % any input parameters you experiments need; here, we have none
            X_data = X{2}{3};  % summary metrics
            
            x_ = obj.rnd([N_samps,M_repet]);  % generate data
            N_2 = N_samps - sum(x_,1);        % number of 0s in the simulated data
            err_2 = abs(N_2 - X_data);

            experiment_pass(:,1) = err_1 <= errors(1,1); % is simulated summary metric sufficiently close?
            experiment_pass(:,2) = err_2 <= errors(1,2); % is simulated summary metric sufficiently close?               
            all_pass = prod(experiment_pass,2);
            approx_L1 = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
            
            experiment_pass(:,1) = err_1 <= errors(2,1,1); % is simulated summary metric sufficiently close?
            experiment_pass(:,2) = err_2 <= errors(2,1,2); % is simulated summary metric sufficiently close?               
            all_pass = prod(experiment_pass,2);
            approx_L2 = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
            
            if approx_L2 == 0
                d_approx_L = 0;
            else
                d_approx_L = approx_L2 / approx_L1;
            end

        end


        

    end





end