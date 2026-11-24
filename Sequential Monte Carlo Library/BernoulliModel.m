%% Bernoulli Model class for testing SMC and ABC SMC classes
%{

%}

classdef BernoulliModel < handle
    %{

    %}

    properties
        p       % probability of a 1
    end

    methods
        function obj = BernoulliModel()
            obj.p = 0.5; % initializing
        end

        function setParams(obj,theta)
            % updating the parameter
            obj.p = theta(1);
        end

        function logL = log_likelihood(obj,X)
            % log likelihood of the dataset X
            N1 = sum(X);    % number of datapoints equal to 1
            N0 = length(X) - N1;    % number of datapoints equal to 0
            
            logL = N1*log(obj.p) + N0*log(1-obj.p);
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
    end
end