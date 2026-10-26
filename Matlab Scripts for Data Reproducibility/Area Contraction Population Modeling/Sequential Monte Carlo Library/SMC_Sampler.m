%% SMC Sampling Class
%{
    Last Updated: MJB 10/25/2026
%}


classdef SMC_Sampler < Dummy_Sampler
%{
    This class performs Sequential Monte Carlo (SMC) parameter estimation 
    for a provided statistical model. This model must be a class instance 
    that contains at least the functions:
        - setParams(theta) which updates the model's parameters using the
        vector theta
        - log_likelihood(X) which returns the log likelihood for a
        dataset X (where each row is a sample of theta)
   
    Note: this implementation of the SMC algorithm only allows for scalar
    parameters, each with distinct proposal kernels (no joint proposal 
    kernels). If your statistical model required vector or matrix 
    parameters, they must be constructed from scalar parameters inside the 
    setParams function.
%}

    properties
        
        gammas          % column vector of exponents used to interpolate between the prior and the posterior
                       
        M_generations   % number of generations over which to sample the parameters (length of gamma)
        propose_mult    % value used to multiple the emperical covarinace when determining the proposal width

        model           % statistical model for which you are estimating the parameters
        proposals       % vector of proposal widths, updated based on particle variance
        
        transforms      % transformer classes that allows for bounded variables to be transformed to the full real number line.
                        % It also calculates the necessary corrections to the prior probability in the MH step

        data            % dataset used to estimate your parameters; each row is assumed to be an independent sample from your population model

        
        particle_history    % matrix of particles used at each generation after resampling; used to show the progression of the sampler
        weight_history      % matrix of particle weights
    end

    methods
        function obj = SMC_Sampler(N,model,priors, propose_mult,data,gammas,parameter_labels,maps,transforms,report_prior)
            arguments
                N 
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
            
            % calculating values to set up the Dummy Sampler
            if isempty(maps)
                P_parameters = length(priors);
            else
                P_parameters = 0;
                for i=1:size(maps,1)
                    P_parameters = P_parameters + length(maps{i,2});
                end
            end
            
            % running super initial function; sets N_particles,
            % P_parameters, maps, priors, report_priors
            posterior_ = zeros(N,P_parameters);
            obj@Dummy_Sampler(posterior_,parameter_labels,maps,priors,report_prior)
           
            % additional parameters
            obj.gammas = gammas;
            obj.M_generations = length(gammas);

            obj.model = model;
            obj.transforms = transforms;

            obj.data = data;
            
            obj.proposals = diag(ones(obj.P_parameters,1));
            obj.propose_mult = propose_mult;

            obj.particle_history = zeros(obj.M_generations,obj.N_particles,obj.P_parameters);   
            obj.weight_history = zeros(obj.M_generations,obj.N_particles);
            obj.calibrated = false;
        end


        % ~~~~~ functions for calibrating the model/estimating the posterior ~~~~~
        function GeneratePosterior(obj)
            %{
                This function runs the full SMC algorithm
            %}  
            
            % preallocating memory for the vectors/matrices used at each
            % stage
            cur_particles = zeros(obj.N_particles,obj.P_parameters);
            weights = ones(obj.N_particles,1) / obj.N_particles;

            % initializing the parameters using the priors
            if isempty(obj.maps)
                for i=1:obj.P_parameters
                    cur_particles(:,i) = obj.priors{i}.random([obj.N_particles,1]);
                end
            else
                for i=1:size(obj.maps,1)
                    prior_ind = obj.maps{i,1};
                    prior_samps = obj.priors{prior_ind}.random([obj.N_particles,1]);
                    particle_inds = obj.maps{i,2};
                    for j=1:length(particle_inds)
                        cur_particles(:,particle_inds(j)) = prior_samps(:,j);
                    end
                end
            end
            
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
                
                % calculate the new weights for particles
                tic;
                weights_updates = obj.CalculateWeightUpdates(cur_particles,kk);
                weights = weights .* weights_updates;   % updating particles;  
                weights = weights / sum(weights);       % normalizing weights
                t1 = toc;

                % calculating the effective sample size (ESS) based on the
                % new weights
                ESS = 1./sum( weights.^2 );
                fprintf("\tESS: %.3f, Full Sample Size: %i, ESS%%: %.3f\n",ESS, obj.N_particles, ESS/obj.N_particles);

                % Effective sample size drops below threshold. Time to
                % resample and mutate.
                resample = (ESS < 0.75*obj.N_particles);
                if resample
                    fprintf("\tESS below threshold (alpha*N_particles)! Resampling...\n")

                    % resample particles based on weights to prevent degeneracy
                    particle_index = obj.SystematicResample(weights,obj.N_particles);
                    cur_particles = cur_particles(particle_index,:);
                    weights = ones(obj.N_particles,1) / obj.N_particles; % resetting weights because particles were weighted resampled
                end

                % update proposal widths based on jointed weighted covariance
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
                

                % perturb current particles
                tic;
                [cur_particles,acceptance_rate] = obj.PerturbParticles(cur_particles,kk);
                t2 = toc;

                fprintf("\tMutation Acceptance Rate: %.1f%%\n",100*acceptance_rate);
                
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

        function [particles,acceptance_rate] = PerturbParticles(obj,cur_particles,generation)
            %{
                Performs an MCMC perturbation step on the current particles
            %}
            particles = 0*cur_particles;    % preallocating
            
            % creating a local function to pass to the parallel loop
            PerturbParticle_loc = @(cur_particle) obj.PerturbParticle_i(cur_particle,obj.proposals,obj.model,obj.data,obj.gammas(generation,:),@(x) obj.CalculatePrior(x),obj.transforms);
            
            function update_progress()
                count = count + 1;
                if mod(count,one_percent)==0
                    % Calculate percentages and bar length
                    bar_length = 20; % Total width of the text bar
                    progress_ticks = round((count / obj.N_particles) * bar_length);
                    
                    % Construct the text string: e.g., "#####               ] 25%"
                    bar_str = [repmat('#', 1, progress_ticks), repmat(' ', 1, bar_length - progress_ticks)];
                    pct_str = sprintf('%3d%%', round((count / obj.N_particles) * 100));
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
            one_percent = floor(obj.N_particles / 100);
            total_backspaces = 0;
            afterEach(dq, @(msg) update_progress());

            accepted_total = 0;
            parfor j = 1:obj.N_particles
               [particles(j,:),accepted] = PerturbParticle_loc(cur_particles(j,:));
               accepted_total = accepted_total + accepted;
               send(dq,true);
            end             
            acceptance_rate = accepted_total / obj.N_particles;
        end
        
        function [particle_i, accepted] = PerturbParticle_i(obj,cur_particle,proposal_widths,model_i,data_,gammas_,prior,transforms)
            %{
                Performs perturbation for a single particle
            %}
            particle_i = cur_particle;
            accepted = 0;
            N_steps = 30;
            for i=1:N_steps
                cur_particle = particle_i;
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
    
    
                % do a Metropolis update
                
                % calculate the log likelihood for the current particles
                model_i.setParams(cur_particle);
                cur_likelihood =  gammas_ * model_i.log_likelihood(data_) + log(prior(cur_particle));   % log likelihood in the current generation
                
                % calculate the log likelihood for the proposed particles
                model_i.setParams(proposed_particle);
                proposed_likelihood =  gammas_ * model_i.log_likelihood(data_) + log(prior(proposed_particle));
                
                % if you have any transformed particles, calculate the jacobian
                % correction to the likelihood ratios
                correction = 0;
                if ~isempty(transforms)
                    for j=1:length(transforms)
                        ind = transforms{j}{1};
                        correction = correction + log(transforms{j}{2}.jacobian(cur_particle(ind))) - log(transforms{j}{2}.jacobian(proposed_particle(ind)));
                    end
                end
                
                % Metropolis Hastings step
                MH_log_ratio = min([0, proposed_likelihood - cur_likelihood + correction]);
                p_accept = rand(1);
                if exp(MH_log_ratio) > p_accept
                    % accept the particle
                    particle_i = proposed_particle;
                    accepted = accepted + 1;
                else
                    % accept the previous particle
                    particle_i = cur_particle;
                    accepted = accepted + 0;
                end
            end
            accepted = accepted / N_steps;
        end

        function weights = CalculateWeightUpdates(obj,particles,generation)
            %{
                Calculates the weights associated with each particles
            %}

            CalculateWeightUpdate_loc = @(particle) obj.CalculateWeightUpdate_i(particle,obj.model,obj.gammas,obj.data,generation);
            weights = zeros(obj.N_particles,1);
            parfor i=1:obj.N_particles
                weights(i) = CalculateWeightUpdate_loc(particles(i,:));
            end
        end

        function weight_i = CalculateWeightUpdate_i(obj,particle_i,model_i,gammas,data,generation)
            %{
                Calculates the weight update for a single particle
            %}

            model_i.setParams(particle_i);
            cur_likelihood = model_i.log_likelihood(data);
            log_weight_i = (gammas(generation) - gammas(generation-1))*cur_likelihood;
            weight_i = exp(log_weight_i);
        end

        

        
        %%%%%% Fixed Helper Functions Used By SMC_Sampler and ABC_SMC_Sampler %%%%%%
        function weighted_cov = CalculateWeightedCovariance(obj,weights,particles)
            % calculates the weighted emperical covariance of the particles
            % using the weights in weights
            weighted_mean = sum(weights .* particles);  % weighted mean
            d_mu = particles - weighted_mean;           % difference from mean
            weighted_cov = (weights.*d_mu)' * d_mu;     % weighted covariance
        end    
        

        function figs = PlotConvergence(obj)
            %{
                Plots the convergence of the different parameters in the
                model
            %}
            if ~obj.calibrated
                error("Posterior must be calibrated before you can plot the convergence")
            end
            
            % for each of the parameters, plot the convergence of the mean
            % and 95% credible interval over the M generations
            figs = cell(obj.P_parameters,1);

            for i=1:obj.P_parameters
                fi = figure("Position",[100,100,600,300],'Color','w'); hold all
                
                for k=1:obj.M_generations
                    weights_p = obj.weight_history(k,:);
                    resamp_inds = obj.SystematicResample(weights_p,obj.N_particles);
                    particles_p = obj.particle_history(k,resamp_inds,i);
                    
                    mean_ = median(particles_p);
                    cred_int = quantile(particles_p,[0.025,0.975]);
                    errorbar(k,mean_,cred_int(1)-mean_,cred_int(2)-mean_,'vertical','LineWidth',1,'Color','k','Marker','s','MarkerFaceColor','k');
                end
                xlabel("Generation [ ]");
                ylabel(obj.parameter_labels(i));
                set(gca,"FontName","calibri","FontSize",15);
                figs{i} = fi;
            end
    
        end

    end



end

