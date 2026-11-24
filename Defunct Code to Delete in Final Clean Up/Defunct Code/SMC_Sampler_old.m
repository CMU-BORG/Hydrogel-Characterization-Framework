%% SMC Sampling Class
%{
    Last Updated: MJB 10/10/2026
%}


classdef SMC_Sampler < handle
%{
    This class performs Sequential Monte Carlo (SMC) parameter estimation 
    for a provided statistical model. This model must be a class instance 
    that contains at least the functions:
        - setParams(theta) which updates the model's parameters using the
        vector theta
        - log_likelihood(X) which returns the log likelihood for a
        dataset X (where each row is a sample of theta)

    The class also has an optional input 'use_ABC' that enables the sampler
    to use Approximate Bayesian Computation to replace the need for an
    explicit likelihood function. This also requires the model class to
    contain a function:
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
        
        gammas          % NOTE: THIS PARAMETER IS OVERLOADED FOR THE STANDARD AND ABC SAMPLERS
                        %   FOR STANDARD SMC: column vector of exponents used to interpolate between the prior and the posterior
                        %   FOR ABC SMC: matrix of error tolerance where each slice contains an error tolerance for each
                        %                of the summary metrics in the supplied dataset. The first index of the matrix is for the
                        %                generation/level of error tolerance, the second is the each different dataset, and the 
                        %                final is for the different metrics of a given dataset (multiple metrics calculated from 
                        %                the same simulated data).

        N_particles     % number of particles to generate
        M_generations   % number of generations over which to sample the parameters (length of gamma)
        propose_mult    % value used to multiple the emperical covarinace when determining the proposal width

        model           % statistical model for which you are estimating the parameters
        priors          % set of priors for your model parameters, must be a list of distributions of the same length as your list of scalar vectors
        proposals       % vector of proposal widths, updated based on particle variance
        maps            % mappings between priors and the parameters of the the particles. If {}, then every prior is used for one particle
                        % if using maps, the first element of the cell array is the index into the priors list, and the second element is the vector 
                        % of parameter indices that the output of the joint prior point to
        
        data            % dataset used to estimate your parameters; For standard SMC and ABC SMC, the structure of data will be different
                        %   FOR STANDARD SMC: each row is assumed to be an independent sample from your population model
                        %   FOR ABC SMC: data is a cell array where each row is a set of summary metrics and experimental parameters


        P_parameters    % number of parameters in the model (calculated from the number of priors)
        parameter_labels    % array of string labels for each the parameters. Used for labeling plots
        
        posterior       % vector of particles that constitute the estimated prior (N_particlesxP_parameters)
        particle_history    % matrix of particles used at each generation after resampling; used to show the progression of the sampler
        weight_history      % matrix of particle weights
        calibrated      % boolean variable set to true once GeneratePosterior has been run
        
        bw  % bw for the kernel density estimation of the posterior based on Silverman's Rule of Thumb
        report_prior % when the pdf is called, the total prior pdf is reported instead of the estimated posterior. This is used for multistep estimations 

        use_ABC         % optional boolean variable set to False (default) for standard Bayesian approaches, or to True for Approximate Bayesian Computation methods
        S_sims          % number of simulations to run to approximate the likelihood
    end

    methods
        function obj = SMC_Sampler(N,model,priors, propose_mult,data,gammas,parameter_labels,use_ABC,maps,report_prior)
            arguments
                N 
                model 
                priors 
                propose_mult 
                data 
                gammas 
                parameter_labels 
                use_ABC = 0;
                maps = [];
                report_prior = 0;
            end


            obj.N_particles = N;
            obj.gammas = gammas;
            obj.M_generations = length(gammas);

            obj.model = model;
            obj.priors = priors;
            obj.maps = maps;
            obj.data = data;
            if isempty(maps)
                obj.P_parameters = length(priors);
            else
                P_parameters = 0;
                for i=1:size(maps,1)
                    P_parameters = P_parameters + length(maps{i,2});
                end
                obj.P_parameters = P_parameters;
            end
            obj.proposals = diag(ones(obj.P_parameters,1));
            obj.parameter_labels = parameter_labels;
            obj.propose_mult = propose_mult;

            obj.posterior = zeros(obj.N_particles,obj.P_parameters);
            obj.particle_history = zeros(obj.M_generations,obj.N_particles,obj.P_parameters);   
            obj.weight_history = zeros(obj.M_generations,obj.N_particles);
            obj.calibrated = false;

            obj.use_ABC = use_ABC;
            obj.S_sims = 20;

            % Settings for using the sampler result as an input for another
            % sampler
            obj.bw = zeros(obj.P_parameters,1); % used to calculate the bandwidth of the KDE
            obj.report_prior = report_prior; 
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

            % for i=1:obj.P_parameters
            %    obj.proposals(i) = sqrt( (obj.propose_mult^2 / obj.P_parameters)*var(cur_particles(:,i)) ); % update proposal widths based on particle variance
            % end
            weighted_mean = sum(weights .* cur_particles);
            d_mu = cur_particles - weighted_mean;
            weighted_cov = (weights.*d_mu)' * d_mu;
            obj.proposals = (obj.propose_mult^2 / obj.P_parameters )*weighted_cov; 

            obj.particle_history(1,:,:) = cur_particles;
            obj.weight_history(1,:) = weights;

            indices = 1:obj.N_particles;

            for kk = 2:obj.M_generations
                % loop over all of the intermediate generations
                fprintf("Starting Generation %i...\n",kk)
                
                % calculate the new weights for particles
                tic;
                weights_updates = obj.CalculateWeightUpdates(cur_particles,kk);
                if obj.use_ABC
                    weights = weights_updates;
                else
                    weights = weights .* weights_updates;   % updating particles;  
                end
                weights = weights / sum(weights);       % normalizing weights
                t1 = toc;

                ESS = 1./sum( weights.^2 );
                fprintf("\tESS: %.3f, Full Sample Size: %i, ESS%%: %.3f\n",ESS, obj.N_particles, ESS/obj.N_particles);

                % Effective sample size drops below threshold. Time to
                % resample and mutate.
                resample = (ESS < 0.5*obj.N_particles) | (obj.use_ABC);
                if resample
                    if obj.use_ABC
                        fprintf("\tResampling...\n")
                    else
                        fprintf("\tESS below threshold (alpha*N_particles)! Resampling...\n")
                    end
                    % resample particles based on weights to prevent degeneracy
                    % particle_index = randsample(indices,obj.N_particles,true,weights);
                    particle_index = obj.SystematicResample(weights,obj.N_particles);
                    cur_particles = cur_particles(particle_index,:);
                    weights = ones(obj.N_particles,1) / obj.N_particles; % resetting weights because particles were weighted resampled
                end
                
                % update proposal widths based on jointed weighted covariance
                weighted_mean = sum(weights .* cur_particles);
                d_mu = cur_particles - weighted_mean;
                weighted_cov = (weights.*d_mu)' * d_mu;
                obj.proposals = (obj.propose_mult^2 / obj.P_parameters )*weighted_cov; 
                
                tic;
                % perturb current particles
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
            
            % creating local copies of variables to passing into the for
            % loop
            model_i = obj.model;
            data_ = obj.data;
            proposal_widths = obj.proposals;
            gammas_ = obj.gammas(generation,:);
            weights = obj.weight_history(generation-1,:);

            % priors_ = obj.priors;

            use_ABC_ = obj.use_ABC;
            
            PerturbParticle_loc = @(cur_particles) obj.PerturbParticle_i(cur_particles,proposal_widths,model_i,data_,gammas_,weights);
            
            function update_progress()
                count = count + 1;
                
                % Calculate percentages and bar length
                bar_length = 20; % Total width of the text bar
                progress_ticks = round((count / obj.N_particles) * bar_length);
                
                % Construct the text string: e.g., "#####               ] 25%"
                bar_str = [repmat('#', 1, progress_ticks), repmat(' ', 1, bar_length - progress_ticks)];
                pct_str = sprintf('] %3d%%', round((count / obj.N_particles) * 100));
                
                % Clear the previous line in the terminal using backspaces (\b)
                % We backspace the exact length of the bar + closing bracket + percentage text
                total_backspaces = bar_length + length(pct_str);
                fprintf(repmat('\b', 1, total_backspaces));
                
                % Print the updated text
                fprintf('%s%s', bar_str, pct_str);
            end

            % 1. Create a DataQueue for worker-to-client communication
            dq = parallel.pool.DataQueue();
            
            % 2. Define a persistent counter to track completed iterations
            % and attach the update function as a listener to the queue
            count = 0;
            afterEach(dq, @(msg) update_progress());


            accepted_total = 0;
            parfor j = 1:obj.N_particles
               if use_ABC_
                   accepted = -1;
                   while accepted < 0
                       % loop until at least one simulation was accepted 
                       [particle_i,accepted] = PerturbParticle_loc(cur_particles);
                   end
               else
                    [particle_i,accepted] = obj.PerturbParticle_i(cur_particles(j,:),proposal_widths,model_i,data_,gammas_,weights);
               end
               particles(j,:) = particle_i; 
               accepted_total = accepted_total + accepted;
               % fprintf("Number Accepted: %i\n",accepted_total);
               send(dq,true);
            end             
            acceptance_rate = accepted_total / obj.N_particles;
        end
        
        function [particle_i, accepted] = PerturbParticle_i(obj,cur_particle,proposal_widths,model_i,data_,gammas_,weights)
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
            % if isempty(obj.maps)
            %     for i = 1:length(proposed_particle)
            % 
            %         cur_i = cur_particle(i);
            % 
            %         while ~valid_point
            %             prop_i = cur_i + proposal_widths(i)*normrnd(0,1,1); 
            %             if priors_{i}.pdf(prop_i) > 0
            %                 % point is valid if the prior probability is
            %                 % greater than 0
            %                 valid_point = 1;
            %             end
            %         end
            %         proposed_particle(i) = prop_i;
            %     end
            % else
            %     for i=1:size(obj.maps,1)
            %         prior_ind = obj.maps{i,1};
            %         particle_inds = obj.maps{i,2};
            % 
            %         valid_point = 0;
            %         cur_i = cur_particle(particle_inds);
            %         prop_i = 0*cur_i;
            %         while ~valid_point
            %             % perturb each particle in the joint prior
            %             for j=1:length(particle_inds)
            %                 prop_i(j) = cur_i(j) + proposal_widths(particle_inds(j))*normrnd(0,1,1); 
            %             end
            %             if priors_{prior_ind}.pdf(prop_i) > 0
            %                 valid_point = 1;
            %             % else
            %             %     fprintf("%s Out of prior range!\n",obj.parameter_labels{prior_ind})
            %             end
            %         end
            %         proposed_particle(particle_inds) = prop_i;
            %     end
            % end

            if obj.use_ABC
                % for every perturbation in ABC, sample a new value
                ind = obj.SystematicResample(weights,1);
                cur_particle = cur_particle(ind,:);
            end

            valid_point = 0;
            while ~valid_point
                if isscalar(cur_particle)
                    proposed_particle = cur_particle + normrnd(0,sqrt(proposal_widths),1);
                else
                    proposed_particle = cur_particle + mvnrnd(0*cur_particle, proposal_widths,1);
                end
                valid_point = obj.CalculatePrior(proposed_particle) > 0;
            end

            if obj.use_ABC
                model_i.setParams(proposed_particle);
                proposed_likelihood = log(model_i.approx_likelihood(data_,gammas_,obj.S_sims));
                if proposed_likelihood == -Inf
                    % no simulations produced sufficiently close data
                    particle_i = cur_particle;
                    accepted = -1;
                else
                    particle_i = proposed_particle;
                    accepted = 1;
                end        
            else
                % do a Metropolis update
                
                % calculate the log likelihood for the current particles
                model_i.setParams(cur_particle);
                cur_likelihood =  gammas_ * model_i.log_likelihood(data_) + log(obj.CalculatePrior(cur_particle));   % log likelihood in the current generation
                
                % calculate the log likelihood for the proposed particles
                model_i.setParams(proposed_particle);
                proposed_likelihood =  gammas_ * model_i.log_likelihood(data_) + log(obj.CalculatePrior(proposed_particle));
                
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
        end

        function weights = CalculateWeightUpdates(obj,particles,generation)
            %{
                Calculates the weights associated with each particles
            %}
            old_weights = obj.weight_history(generation-1,:);
            old_particles = squeeze(obj.particle_history(generation-1,:,:));
            weights = zeros(obj.N_particles,1);
            
            % creating local copies of variables to passing into the for
            % loop
            model_i = obj.model;
            gammas_ = obj.gammas;
            data_ = obj.data;
            proposal_widths_ = obj.proposals;
            parfor i=1:obj.N_particles
                weights(i) = obj.CalculateWeightUpdate_i(particles(i,:),model_i,gammas_,data_,old_particles,generation,old_weights,proposal_widths_);
            end
        end

        function weight_i = CalculateWeightUpdate_i(obj,particle_i,model_i,gammas,data,old_particles,generation,old_weights,proposal_widths)
            %{
                Calculates the weight update for a single particle
            %}

            model_i.setParams(particle_i);
            if obj.use_ABC
                % cur_likelihood = model_i.approx_likelihood(data,gammas(generation,:)) ;  % tolerance added to prevent divide by 0
                % prev_likelihood = model_i.approx_likelihood(data,gammas(generation-1,:)) ;
                % if cur_likelihood == 0
                %     weight_i = 0;
                % else
                %     weight_i = cur_likelihood / prev_likelihood; 
                % end

                % gammas_ = zeros([2,size(gammas(generation,:))]);
                % gammas_(1,:,:) = gammas(generation-1,:);
                % gammas_(2,:,:) = gammas(generation,:);
                
                % get the mutliplier based on how many simulations lead to
                % sufficiently close data if this set of simulations led to 
                % 0 weight, we reset it to 1/S_sims, because we know it can 
                % generate at least 1 successful simulation from the perturbation
                b_i = max(model_i.approx_likelihood(data, gammas(generation,:),obj.S_sims), 1/obj.S_sims);
                if generation == 1
                    weight_i = b_i;
                else
                    % d_particle = particle_i - old_particles;
                    if isscalar(particle_i)
                        log_trans = log(normpdf(old_particles,particle_i,sqrt(proposal_widths))');
                    else
                        log_trans = log(mvnpdf(old_particles,particle_i,proposal_widths));
                    end
                    log_transition_prob = log(old_weights)' + log_trans;
                    % for i=1:obj.P_parameters
                    %     log_transition_prob = log_transition_prob + log(normpdf(particle(i) - particles(:,i), 0, proposal_widths(i)))';
                    % end
                    weight_i = obj.CalculatePrior(particle_i) * b_i / sum( exp(log_transition_prob) );
                end
            else
                cur_likelihood = model_i.log_likelihood(data);
                log_weight_i = (gammas(generation) - gammas(generation-1))*cur_likelihood;
                weight_i = exp(log_weight_i);
            end
            end

        function p_ = CalculatePrior(obj,particles)
            % calculates the total prior probability of a particle
            if isempty(obj.maps)
                log_pdf_ = zeros(size(particles));
                for i=1:obj.P_parameters
                    log_pdf_(:,i) = log(obj.priors{i}.pdf(particles(:,i)));
                end
            else
                log_pdf_ = zeros(size(particles,1),size(obj.maps,1));
                for i=1:size(obj.maps,1)
                    prior_ind = obj.maps{i,1};
                    particle_inds = obj.maps{i,2};
                    mapped_particles = particles(:,particle_inds);
                    log_pdf_(:,i) = log(obj.priors{prior_ind}.pdf(mapped_particles));
                end
            end
            p_ = exp( sum(log_pdf_,2) );
        end

        function inds = SystematicResample(obj,weights,M)
            %{
                Samples M indices using Systematic Resampling where the
                indices are weighted by weights
            %}

            W_cumul = cumsum(weights);
            U = rand(1)/M;  % offset
            u_i = U + (0:M-1)/M;
            inds = zeros(M,1);
            cur_segment = 1;
            for i=1:M
                while u_i(i)>W_cumul(cur_segment)
                    cur_segment = cur_segment + 1;
                end
                inds(i) = cur_segment;
            end


        end
        
        % ~~~~ functions for generating samples from the posterior ~~~~
        function samps = random(obj,size)
            %{
                Generating random samples from the posterior. Allows for
                samplers to be used for priors of other samplers in
                multistage estimation routines
            %}      
            len = prod(size);
            % indices = 1:obj.N_particles;
            % particle_index = randsample(indices,len,true);
            particle_index = obj.SystematicResample(ones(obj.N_particles,1)/obj.N_particles,len);
            cur_particles = obj.posterior(particle_index,:);
            if obj.report_prior
                samps = cur_particles; % don't perturb the particles if you are using this sampler as an intermediate distribution
            else
                % if this is the finalized distribution, adding this
                % perturbation step may help mixing
                samps = obj.PerturbParticles(cur_particles,obj.M_generations);
            end
        end

        function p_ = pdf(obj,particles)
            %{
                Estimate the pdf using posterior particles and kernel
                density estimation. Uses priors to ensure that disallowed
                particles recieve 0 pdf.
            %}
            p_prior = obj.CalculatePrior(particles);
            if obj.report_prior
                % report the total prior probability for multistep
                % estimation
                p_ = p_prior;
            else
                % estimate the posterior pdf using kernel density
                % estimation
                p_ = mvksdensity(obj.posterior,particles,'Bandwidth',obj.bw);
                    
                is_range = p_prior > 0;
                p_ = p_ .* is_range;
            end
        end
        
        % ~~~~~ functions for analyzing the posterior and convergence ~~~~~
        function summary = PosteriorSummary(obj)
            %{
                Creates a struct with the summary metrics for your
                posterior
            %}  
            if ~obj.calibrated
                error("Posterior must be calibrated before you can plot the convergence")
            end
            
            summary = zeros(obj.P_parameters,3);
            fprintf("~~~~ Posterior Summary (mean [95%% CI]) ~~~~\n")
            for p=1:obj.P_parameters
                summary(p,1) = mean(obj.posterior(:,p));
                summary(p,2:3) = quantile(obj.posterior(:,p),[0.025,0.975]);
                fprintf("\t%s: %.3f [ %.3f,%.3f ] \n",obj.parameter_labels{p},summary(p,1),summary(p,2),summary(p,3))
            end

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
            indices = 1:obj.N_particles;
            for i=1:obj.P_parameters
                fi = figure("Position",[100,100,600,300],'Color','w'); hold all
                
                for k=1:obj.M_generations
                    weights_p = obj.weight_history(k,:);
                    % resamp_inds = randsample(indices,obj.N_particles,true,weights_p); % resampling to compensate for weighting
                    resamp_inds = obj.SystematicResample(weights_p,obj.N_particles);
                    particles_p = obj.particle_history(k,resamp_inds,i);
                    
                    mean_ = mean(particles_p);
                    cred_int = quantile(particles_p,[0.025,0.975]);
                    errorbar(k,mean_,cred_int(1)-mean_,cred_int(2)-mean_,'vertical','LineWidth',1,'Color','k','Marker','s','MarkerFaceColor','k');
                end
                xlabel("Generation [ ]");
                ylabel(obj.parameter_labels(i));
                set(gca,"FontName","calibri","FontSize",15);
                figs{i} = fi;
            end
    
        end

        function figs = PlotPosteriorHistograms(obj)
            %{
                Plots the histograms of particles for each parameter
            %}
            if ~obj.calibrated
                error("Posterior must be calibrated before you can plot the convergence")
            end
            
            % for each of the parameters, plot the histograms of particles
            figs = cell(obj.P_parameters,1);
            for i=1:obj.P_parameters
                fi = figure("Position",[100,100,600,300],'Color','w'); hold all
                particles_p = obj.posterior(:,i); 
                histogram(particles_p,'Normalization','pdf')
                xlabel(obj.parameter_labels(i));
                set(gca,"FontName","calibri","FontSize",12);
                figs{i} = fi;
            end

        end

        function fig = PlotJointPosterior(obj)
            fig = figure("Position",[100,100,850,850],'Color','w');
            tiledlayout(obj.P_parameters,obj.P_parameters,"TileSpacing","tight","Padding","compact"); 
            
            for i=1:obj.P_parameters
                ind = sub2ind([obj.P_parameters,obj.P_parameters],i,i);
                axi = nexttile(ind);
                histogram(axi,obj.posterior(:,i),"Normalization","pdf");
                title(obj.parameter_labels{i})
                
                for j=i+1:obj.P_parameters
                    ind = sub2ind([obj.P_parameters,obj.P_parameters],i,j);
                    axi = nexttile(ind);
                    plot(axi, obj.posterior(:,i),obj.posterior(:,j),'.k')
                end
                    
            end
            

        end

    end



end

