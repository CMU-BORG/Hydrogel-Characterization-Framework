%% ABC-SMC Sampling Class
%{
    Last Updated: MJB 10/7/2026
%}

classdef ABC_SMC_Sampler < SMC_Sampler
%{
    This class modified the Sequential Monte Carlo (SMC) Sampler class to
    utilize Approximate Bayesian Computation instead of traditional
    likelihood caculations. As with the SMC_Sampler class, this class
    performs parameter estimation for a provided statistical model. This 
    model must be a class instance that contains at least the functions:
        - setParams(theta) which updates the model's parameters using the
        vector theta
        - 

    Note: this implementation of the ABC-SMC algorithm only allows for 
    scalar parameters, each with distinct priors (no joint priors) and 
    distinct proposal kernels (no joint proposal kernels). If your 
    statistical model required vector or matrix parameters, they must be 
    constructed from scalar parameters inside the setParams function.
%}

    properties

        

    end

    methods
        function obj = ABC_SMC_Sampler(N,model,priors,propose_mult,data,gammas,parameter_labels)
            obj@SMC_Sampler(N,model,priors, propose_mult,data,gammas,parameter_labels);
            
        end

        % ~~~~~ functions for calibrating the model/estimating the posterior ~~~~~
        % some functions are identical to the those in the SMC Sampler 
        %   - functions that remain the same:
        %       - GeneratePosterior
        %       - PerturbParticles

        
        function [particle_i, accepted] = PerturbParticle_i(obj,cur_particle,proposal_widths,priors_,model_i,data_)
            %{
                Performs perturbation for a single particle
            %}

            % for each particle, propose a new value based on the
            % perturbation kernel, check if it is within support of the
            % priors, and perform a Metropolis-Hastings check to see if
            % the proposed particleis accepted
            proposed_particle = 0*cur_particle;
            
            % perturb each parameter in the particle so that they are
            % all valid points
            for i = 1:length(proposed_particle)
                valid_point = 0;
                cur_i = cur_particle(i);

                while ~valid_point
                    prop_i = cur_i + proposal_widths(i)*normrnd(0,1,1); %proposals_(i).random(1);
                    if priors_(i).pdf(prop_i) > 0
                        % point is valid if the prior probability is
                        % greater than 0
                        valid_point = 1;
                    end
                end
                proposed_particle(i) = prop_i;
            end

            % calculate the log likelihood for the current particles
            % model_i = copy(obj.model);
            model_i.setParams(cur_particle);
            cur_likelihood = model_i.log_likelihood(data_);
            
            % calculate the log likelihood for the proposed particles
            model_i.setParams(proposed_particle);
            proposed_likelihood = model_i.log_likelihood(data_);
            
            % Metropolis Hastings step
            MH_log_ratio = min([0, proposed_likelihood - cur_likelihood]);
            p_accept = rand(1);
            if exp(MH_log_ratio) > p_accept
                % accept the particle
                particle_i = proposed_particle;
                accepted = 1;
            else
                % accept the previous particle
                particle_i = cur_particle;
                accepted = 0;
            end
        end

        function weights = CalculateWeightUpdates(obj,particles,generation)
            %{
                Calculates the weights associated with each particles
            %}

            weights = zeros(obj.N_particles,1);
            
            % creating local copies of variables to passing into the parfor
            % loop
            model_i = obj.model;
            gammas_ = obj.gammas;
            data_ = obj.data;

            parfor i=1:obj.N_particles
                weights(i) = obj.CalculateWeightUpdate_i(model_i,gammas_,data_,particles(i,:),generation);
            end
        end

        function weight_i = CalculateWeightUpdate_i(obj,model_i,gammas,data,particle,generation)
            %{
                Calculates the weight update for a single particle
            %}
            model_i.setParams(particle);
            cur_likelihood = model_i.log_likelihood(data);
            log_weight_i = (gammas(generation) - gammas(generation-1))*cur_likelihood;
            weight_i = exp(log_weight_i);
        end

    end




end