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
