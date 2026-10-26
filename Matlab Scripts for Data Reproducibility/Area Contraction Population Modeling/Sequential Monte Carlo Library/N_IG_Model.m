%% Normal-Inverse-Gamma Model class for testing SMC and ABC SMC classes
%{

%}

classdef N_IG_Model < handle
    %{

    %}

    properties
        mu0       % mean
        v0        % 
        alpha    %
        beta     %
    end

    methods
        function obj = N_IG_Model()
            obj.mu0 = 0; % initializing
            obj.v0 = 1;
            obj.alpha = 1;
            obj.beta = 1;
        end

        function setParams(obj,theta)
            % updating the parameter
            obj.mu0 = theta(1);
            obj.v0 = theta(2);
            obj.alpha = theta(3);
            obj.beta = theta(4);
        end

        function logL = log_likelihood(obj,X)
            % log likelihood of the dataset X
            log_p = log(obj.pdf(X));
            logL = sum(log_p);
        end

        function p = pdf(obj,X)
            mu = X(:,1);
            sigma2 = X(:,2);

            p_sig = gampdf(1./sigma2, obj.alpha,1/obj.beta) .* (1 ./ sigma2.^2);
            p_mu = 0*p_sig;
            for i=1:length(mu)
                p_mu(i) = normpdf(mu(i),obj.mu0,sqrt(sigma2(i)/obj.v0));
            end
            p = p_sig.*p_mu;
        end

        function X = random(obj,N)
            % generates N random Bernoulli variables
            if (length(N)==1)||(N(2)==1)
                sigma2 = 1 ./ gamrnd(obj.alpha, 1/obj.beta, [N(1), 1]); % sampling variance from an inverse gamma
                mu = 0*sigma2;
                for i=1:N(1)
                    mu(i) = normrnd(obj.mu0,sqrt(sigma2(i)/obj.v0),1);
                end
                X = [mu, sigma2];
            else
                error("Multidimensional sampling of this distribution is not available.")
            end
        end
    end
end

        % function approx_L = approx_likelihood(obj,X,errors)
        %     % an approximate likelihood for this model to test the ABC SMC
        %     % code. In this approximate method, N samples are generated
        %     % from the model, and the total numbers of 1s are calculated as
        %     % the summary metric. Then the percent error between the
        %     % simulation and experimental summary metrics is calculated and
        %     % compared to the error threshold. This is conducted for 20
        %     % repetitions
        % 
        %     M_repet = 300; % number of experiments to run; this can be tuned
        % 
        %     % From one dataset get the mean and standard deviation
        %     N_samps = X{1}{1}; % sample size of the experimental data 
        %     Sim_params = X{1}{2};  % any input parameters you experiments need; here, we have none
        %     X_data = X{1}{3};  % summary metrics
        % 
        %     x_ = obj.rnd([N_samps,M_repet]);  % generate all data
        %     metrics = [mean(x_,1)', std(x_,[],1)'];   % calculate the metrics for each of the repeats
        %     errs = (metrics - X_data)./X_data;
        % 
        %     experiment_pass = abs(errs) < errors(1,:);
        %     all_pass = prod(experiment_pass,2);
        % 
        %     approx_L = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
        % end
        % 
        % function d_approx_L = delta_approx_likelihood(obj,X,errors)
        %     % Here, the change in the approximate likelihood between two
        %     % different generations is calculated. Used in weight updates
        % 
        %     M_repet = 300; % number of experiments to run; this can be tuned
        % 
        %     % From one dataset get the mean and standard deviation
        %     N_samps = X{1}{1}; % sample size of the experimental data 
        %     Sim_params = X{1}{2};  % any input parameters you experiments need; here, we have none
        %     X_data = X{1}{3};  % summary metrics
        % 
        %     x_ = obj.rnd([N_samps,M_repet]);  % generate all data
        %     metrics = [mean(x_,1)', std(x_,[],1)'];   % calculate the metrics for each of the repeats
        %     errs = (metrics - X_data)./X_data;
        % 
        %     experiment_pass = abs(errs) < errors(1,:);
        %     all_pass = prod(experiment_pass,2);
        %     approx_L1 = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
        % 
        %     experiment_pass = abs(errs) < errors(2,:);
        %     all_pass = prod(experiment_pass,2);
        %     approx_L2 = sum(all_pass) / M_repet;  % approx likelihood = what fraction of simulations were sufficiently close enough
        % 
        %     if approx_L2 == 0
        %         d_approx_L = 0;
        %     else
        %         d_approx_L = approx_L2 / approx_L1;
        %     end
        % end
%     end
% 
% end