%% ABC SMC Sampling Class
%{
    Last Updated: MJB 10/25/2026
%}


classdef ABC_SMC_Sampler < SMC_Sampler
%{
    This class performs Approximate Bayesian Computation Sequential Monte 
    Carlo (SMC) parameter estimation for a provided statistical model. This 
    model must be a class instance that contains at least the functions:
        - setParams(theta) which updates the model's parameters using the
        vector theta
        - approx_likelihood(X,err) which returns the approximate likelihood for
        a dataset X of summary metrics (cell array) where each row contains 
            - Cell 1: The sample size of the dataset use to calculate the 
              summary metrics
            - Cell 2: a set of simulation parameters needed to reproduce
            the experiment summarized in the summary metric
            - Cell 3: a set of summary metrics from the experiment to
            compare with simulation data
    and err is the vector of error tolerances for each of the summary
    metrics.
   
    Note: this implementation of the SMC algorithm only allows for scalar
    parameters, each with distinct proposal kernels (no joint proposal 
    kernels). If your statistical model required vector or matrix 
    parameters, they must be constructed from scalar parameters inside the 
    setParams function.
%}

    properties
        S_sims  % number of replicates to generate each time you test a new particle


    end

    methods
        function obj = ABC_SMC_Sampler(N,S_sims,model,priors,propose_mult,data,gammas,parameter_labels,maps,transforms,report_prior)
            arguments
                N 
                S_sims
                model 
                priors 
                propose_mult 
                data 
                gammas 
                parameter_labels 
                maps = {};
                transforms = {};
                report_prior = 0;
            end
            obj@SMC_Sampler(N,model,priors,propose_mult,data,gammas,parameter_labels,maps,transforms,report_prior);

            obj.S_sims = S_sims;

        end

        % ~~~~~ functions for calibrating the model/estimating the posterior ~~~~~
        function GeneratePosterior(obj)
            %{
                This function runs the full ABC SMC algorithm
            %} 
            
            % preallocating memory for the vectors/matrices used at each
            % stage
            weights = ones(obj.N_particles,1) / obj.N_particles;

            % initializing the parameters using the priors
            cur_particles = obj.SamplePrior(obj.N_particles);
            
            % calculating the covariance of the proposal matrix
            if isempty(obj.transforms)
                weighted_cov = obj.CalculateWeightedCovariance(weights,cur_particles);
            else
                transformed_particles = cur_particles;
                for j=1:length(obj.transforms)
                    ind = obj.transforms{j}{1};
                    transformed_particles(:,ind) = obj.transforms{j}{2}.transform(cur_particles(:,ind));
                end
                weighted_cov = obj.CalculateWeightedCovariance(weights,transformed_particles);                
            end
            obj.proposals = (obj.propose_mult^2 / obj.P_parameters )*weighted_cov; 

            obj.particle_history(1,:,:) = cur_particles;
            obj.weight_history(1,:) = weights;
            
            for kk = 2:obj.M_generations
                % loop over all of the intermediate generations
                fprintf("Starting Generation %i...\n",kk)
                
                % perturb the particles and calclate b_i
                tic
                past_particles = cur_particles; % saving these for use in the sampler and weight updates
                [cur_particles,bs,generated_samples] = PerturbParticles(obj, weights, past_particles, kk);
                t2 = toc;
                fprintf("\tNumber of Generated Samples: %i, (%% of Generated Samples Accepted: %.2f %%)\n",generated_samples,100*obj.N_particles/generated_samples)
                
                % calculate the new weights of the generated particles
                tic
                weights = obj.CalculateWeightUpdates(cur_particles, weights, past_particles, bs);
                t1 = toc;

                % calculating the effective sample size (ESS) based on the
                % new weights
                ESS = 1./sum( weights.^2 );
                fprintf("\tESS: %.3f, Full Sample Size: %i, ESS%%: %.3f\n",ESS, obj.N_particles, ESS/obj.N_particles);

                % calculating the new weighted covariance matrix
                if isempty(obj.transforms)
                    weighted_cov = obj.CalculateWeightedCovariance(weights,cur_particles);
                else
                    transformed_particles = cur_particles;
                    for j=1:length(obj.transforms)
                        ind = obj.transforms{j}{1};
                        transformed_particles(:,ind) = obj.transforms{j}{2}.transform(cur_particles(:,ind));
                    end
                    weighted_cov = obj.CalculateWeightedCovariance(weights,transformed_particles);                
                end
                obj.proposals = (obj.propose_mult^2 / obj.P_parameters )*weighted_cov; 

                fprintf("\tTime to calculate Weights: %.3f s\n",t1)
                fprintf("\tTime to Perturb Particles: %.3f s\n",t2)
                fprintf("\tApprox Time per Generation: %.3f s\n",t1+t2)
               
                % save particles in the history
                obj.particle_history(kk,:,:) = cur_particles;
                obj.weight_history(kk,:) = weights;
            end
            
            % perform one final resampling to get an unweighted estimate of
            % the posterior
            % particle_index = randsample(indices,obj.N_particles,true,weights);
            particle_index = obj.SystematicResample(weights,obj.N_particles);
            cur_particles = cur_particles(particle_index,:);
            obj.posterior = cur_particles; % same the last set of particles as the posterior  

            bw_scale = (4 / ((obj.P_parameters + 2)*obj.N_particles) )^( 1/(obj.P_parameters + 4) );
            for i=1:obj.P_parameters
                obj.bw(i) = std(obj.posterior(:,i))*bw_scale;
            end
            
            obj.calibrated = true;
        end

        function [particles,bs,generated_samples] = PerturbParticles(obj,weights,cur_particles,generation)
            %{
                This function perturbs particles and calculates simulated
                data using those perturbed particles. If the simulated data
                is sufficiently close to the real data, the b_i is
                incremented. If b_i is greater than 0, then the new
                particle is accepted. Otherwise, we resample from the
                previous generation and try again until b_i > 0
            %}
            
            % the input to this function isn't actually important because
            % it always starts with resampling a particle from the previous
            % generation
            PerturbParticle_loc = @(i) obj.PerturbParticle_i(obj.S_sims, obj.data, cur_particles, obj.model, obj.gammas(generation,:), obj.proposals, obj.transforms, @(n) obj.SystematicResample(weights, n));
            particles = 0*cur_particles;
            bs = zeros(obj.N_particles,1);
            generated_samples = 0;

            
            function update_progress()
                count = count + 1;
                if mod(count,one_percent)==0
                    % Calculate percentages and bar length
                    bar_length = 50; % Total width of the text bar
                    progress_ticks = round((count / obj.N_particles) * bar_length);
                    
                    % Construct the text string: e.g., "#####               ] 25%"
                    bar_str = [repmat('#', 1, progress_ticks), repmat(' ', 1, bar_length - progress_ticks)];
                    pct_str = sprintf('%.1f%%', round((count / obj.N_particles) * 100,1));
                    new_bar = sprintf('\tParticles Perturbed: [%s] %s\n', bar_str, pct_str);
                    % Clear the previous line in the terminal using backspaces (\b)
                    % We backspace the exact length of the bar + closing bracket + percentage text
                    fprintf(repmat('\b', 1, total_backspaces));
                    
                    % Print the updated text
                    fprintf('%s',new_bar);
                    total_backspaces = length(new_bar);
                end
            end

            % 1. Create a DataQueue for worker-to-client communication
            dq = parallel.pool.DataQueue();
            
            % 2. Define a persistent counter to track completed iterations
            % and attach the update function as a listener to the queue
            count = 0;
            one_percent = 0.1*floor(obj.N_particles / 100);
            total_backspaces = 0;
            afterEach(dq, @(msg) update_progress());

            parfor i=1:obj.N_particles
                [particles(i,:), bs(i), generated_samples_i] = PerturbParticle_loc(0);
                generated_samples = generated_samples + generated_samples_i;
                send(dq,true);
            end
        end

        function [particle_i,b_i,generated_samples_i] = PerturbParticle_i(obj, S_sims, data_, cur_particles, model_i, errors, proposal_widths, transforms, resample_func)
            %{
                Performs perturbation for a single particles
            %}

            accepted = 0; 
            generated_samples_i = 0;
            while ~accepted
                % repeat choosing perturbed particles until one generates
                % sufficient samples
                generated_samples_i = generated_samples_i + 1;

                proposed_index = resample_func(1);
                cur_particle = cur_particles(proposed_index, :); 

                if isscalar(cur_particle)
                    proposal_step = normrnd(0,sqrt(proposal_widths),1);
                else
                    proposal_step = mvnrnd(0*cur_particle, proposal_widths,1);
                end
                
                proposed_particle = cur_particle + proposal_step;
                if ~isempty(transforms)
                    % correct the values for the transformed particles
                    for j=1:length(transforms)
                        ind = transforms{j}{1};
                        transformed_proposal = transforms{j}{2}.transform(cur_particle(ind)) + proposal_step(ind); % transform the particle and add the transformed perturbation
                        proposed_particle(ind) = transforms{j}{2}.inv_transform(transformed_proposal);
                    end
                end

                % update the model parameters and run simulations
                model_i.setParams(proposed_particle);
                b_ = model_i.approx_likelihood(data_, errors, S_sims);
                accepted = b_ > 0;  % if you got successful simulations, break out
            end
            particle_i = proposed_particle;
            b_i = b_;
        end

        function new_weights = CalculateWeightUpdates(obj,cur_particles, weights, past_particles, bs)
            %{
                Calculates the updated weights for the new particles given
                the number of successful simulations that they generated
                (bs).
            %}
            
            new_weights = bs .* obj.CalculatePrior(cur_particles);
            CalculateWeightUpdate_loc = @(particle) obj.CalculateWeightUpdate_i(particle, weights, past_particles, obj.proposals, obj.transforms);
            parfor i=1:obj.N_particles
                new_weights(i) = new_weights(i) * CalculateWeightUpdate_loc(cur_particles(i,:));
            end
            new_weights = new_weights / sum(new_weights);   % normalize
        end

        function weight_i = CalculateWeightUpdate_i(obj, particle_i, weights, past_particles, proposal_widths, transforms)
            %{
                For a single particle (particle_i) calculate the weight
                update that multiplies b_i and the prior to get the 
                unnormalized weight of the new particle
            %}  
            
            corrections = 0; % possible corrections to the transition probability cause by transformation
            if ~isempty(transforms)
                % correct the values for the particle
                for j=1:length(transforms)
                    ind = transforms{j}{1};

                    % calculate the jacobian and then transform
                    corrections = corrections + log( transforms{j}{2}.jacobian(particle_i(ind)) );
                    particle_i(ind) = transforms{j}{2}.transform(particle_i(ind));
                    
                    % transform previous particles
                    past_particles(:,ind) = transforms{j}{2}.transform(past_particles(:,ind));
                end
            end
            % log (perturbation kernel * corrections from jacobian * previous weights)
            if isscalar(particle_i)
                log_K = log(normpdf(past_particles, particle_i, sqrt(proposal_widths) )) + corrections + log(weights);
            else
                log_K = log(mvnpdf(past_particles, particle_i, proposal_widths)) + corrections + log(weights);
            end
            weight_i = 1 / sum(exp(log_K));
        end

    end
end