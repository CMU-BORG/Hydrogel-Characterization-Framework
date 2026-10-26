%% Gaussian Model class for testing SMC and ABC SMC classes
%{

%}

classdef GaussianModel < handle
    %{

    %}

    properties
        mu       % mean
        sig       % standard deviation

        input_variance  % specifies if the second element of theta is the standard deviation (0) or the variance (1, for use with conjugate priors)
    end

    methods
        function obj = GaussianModel(input_variance)
            obj.mu = 0; % initializing
            obj.sig = 1;
            obj.input_variance = input_variance;
        end

        function setParams(obj,theta)
            % updating the parameter
            obj.mu = theta(1);
            if obj.input_variance
                obj.sig = sqrt(theta(2));
            else
                obj.sig = theta(2);
            end
        end

        function logL = log_likelihood(obj,X)
            % log likelihood of the dataset X
            log_p = log(normpdf(X,obj.mu,obj.sig));
            logL = sum(log_p);
        end

        function p = pdf(obj,X)
            p = normpdf(X,obj.mu,obj.sig);
        end

        function X = rnd(obj,N)
            % generates N random Bernoulli variables
            if size(N)==1
                X = normrnd(obj.mu,obj.sig,N,1);
            else
                X = normrnd(obj.mu,obj.sig,N);
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

            % From one dataset get the mean and standard deviation
            N_samps = X{1}{1}; % sample size of the experimental data 
            Sim_params = X{1}{2};  % any input parameters you experiments need; here, we have none
            X_data = X{1}{3};  % summary metrics

            x_ = obj.rnd([N_samps,M_repet]);  % generate all data
            metrics = [mean(x_,1)', std(x_,[],1)'];   % calculate the metrics for each of the repeats
            errs = (metrics - X_data)./X_data;

            experiment_pass = abs(errs) < errors(1,:);
            all_pass = prod(experiment_pass,2);

            approx_L = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
        end
    end

end