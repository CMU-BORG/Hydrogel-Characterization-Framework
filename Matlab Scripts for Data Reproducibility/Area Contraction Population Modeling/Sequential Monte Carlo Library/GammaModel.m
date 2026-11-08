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
            arguments
                MS_param 
            end
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
            % generates N random Gamma variables
            if size(N)==1
                X = gamrnd(obj.alpha,obj.beta,N,1);
            else
                X = gamrnd(obj.alpha,obj.beta,N);
            end
        end

    end

end