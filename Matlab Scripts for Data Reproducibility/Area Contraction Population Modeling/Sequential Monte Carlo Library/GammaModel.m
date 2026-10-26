%% Gamma Model class for testing SMC and ABC SMC classes
%{

%}

classdef GammaModel < handle
    %{

    %}

    properties
        alpha       % shape param
        beta        % scale param
        mean_sig_parameterized % set to 1 if you want to parameterize the Gamma model by the mean and standard deviation instead of a and b
    end

    methods
        function obj = GammaModel(MS_param)
            obj.alpha = 1; % initializing
            obj.beta = 1;
            obj.mean_sig_parameterized = MS_param;
        end

        function setParams(obj,theta)
            % updating the parameter
            if obj.mean_sig_parameterized
                mu = theta(1); sig = theta(2);
                obj.alpha = (mu/sig)^2; 
                obj.beta = sig^2 / mu;
            else
                obj.alpha = theta(1);
                obj.beta = theta(2);
            end
        end

        function logL = log_likelihood(obj,X)
            % log likelihood of the dataset X
            log_p = log(gampdf(X,obj.alpha,obj.beta));
            logL = sum(log_p);
        end

        function p = pdf(obj,X)
            p = gampdf(X,obj.alpha,obj.beta);
        end

        function X = rnd(obj,N)
            % generates N random Bernoulli variables
            if size(N)==1
                X = gamrnd(obj.alpha,obj.beta,N,1);
            else
                X = gamrnd(obj.alpha,obj.beta,N);
            end
        end

        function approx_L = approx_likelihood(obj,X,errors,M_repet)
            % an approximate likelihood for this model to test the ABC SMC
            % code. In this approximate method, N samples are generated
            % from the model, and the total numbers of 1s are calculated as
            % the summary metric. Then the percent error between the
            % simulation and experimental summary metrics is calculated and
            % compared to the error threshold. This is conducted for 20
            % repetitions
            
            % M_repet = 300; % number of experiments to run; this can be tuned
            
            % From one dataset get the mean and standard deviation
            N_samps = X{1}{1}; % sample size of the experimental data 
            Sim_params = X{1}{2};  % any input parameters you experiments need; here, we have none
            X_data = X{1}{3};  % summary metrics
            
            x_ = obj.rnd([N_samps,M_repet]);  % generate all data
            metrics = [prod(x_,1)', sum(x_,1)'];   % calculate the metrics for each of the repeats
            errs = (metrics - X_data)./X_data;

            experiment_pass = abs(errs) < errors(1,:);
            all_pass = prod(experiment_pass,2);
    
            approx_L = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
        end

        function d_approx_L = delta_approx_likelihood(obj,X,errors)
            % Here, the change in the approximate likelihood between two
            % different generations is calculated. Used in weight updates

            M_repet = 300; % number of experiments to run; this can be tuned
            
            % From one dataset get the mean and standard deviation
            N_samps = X{1}{1}; % sample size of the experimental data 
            Sim_params = X{1}{2};  % any input parameters you experiments need; here, we have none
            X_data = X{1}{3};  % summary metrics
            
            x_ = obj.rnd([N_samps,M_repet]);  % generate all data
            metrics = [mean(x_,1)', std(x_,[],1)'];   % calculate the metrics for each of the repeats
            errs = (metrics - X_data)./X_data;

            experiment_pass = abs(errs) < errors(1,:);
            all_pass = prod(experiment_pass,2);
            approx_L1 = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
            
            experiment_pass = abs(errs) < errors(2,:);
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