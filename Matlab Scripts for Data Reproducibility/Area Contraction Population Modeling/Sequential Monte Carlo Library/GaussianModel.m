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
    
    end
end