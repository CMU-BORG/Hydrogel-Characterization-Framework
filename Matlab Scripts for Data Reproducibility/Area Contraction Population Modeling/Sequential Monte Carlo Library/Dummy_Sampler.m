%% Dummy Sampler
%{
    Last Updated: MJB 10/25/2026

%}

classdef Dummy_Sampler < handle
%{
    This class serves as a dummy sampler to hold the results of a different
    SMC/ABC-SMC sampler for use as a prior.
%}

    properties
        posterior       % unweighted particles that are used to estimate the sampling distribution
        
        N_particles     % number of particles to generate
        P_parameters    % number of parameters in the model (calculated from the number of priors)

        maps            % mappings between priors and the parameters of the the particles. If {}, then every prior is used for one particle
                        % if using maps, the first element of the cell array is the index into the priors list, and the second element is the vector 
                        % of parameter indices that the output of the joint prior point to

        parameter_labels    % array of string labels for each the parameters. Used for labeling plots

        priors          % set of priors for the parameters, must be a list of distributions of the same length as your list of scalar vectors
        report_prior    % if a set of priors is provided, then instead of estimating the pdf based on the particles, the prior probability is reported
        
        bw  % bw for the kernel density estimation of the posterior based on Silverman's Rule of Thumb
        calibrated      % boolean variable set to true once GeneratePosterior has been run

    end

    methods
        function obj = Dummy_Sampler(posterior,parameter_labels,maps,priors,report_prior)
            arguments 
                posterior 
                parameter_labels
                maps = {};
                priors = {};
                report_prior = 0;
            end

            obj.posterior = posterior;

            obj.parameter_labels = parameter_labels;

            obj.N_particles = size(posterior,1);
            obj.P_parameters = size(posterior,2);

            obj.maps = maps;
            obj.priors = priors;
            obj.report_prior = report_prior & ~isempty(priors);

            bw_scale = (4 / ((obj.P_parameters + 2)*obj.N_particles) )^( 1/(obj.P_parameters + 4) );
            for i=1:obj.P_parameters
                obj.bw(i) = std(obj.posterior(:,i))*bw_scale;
            end
            obj.calibrated = 1;
        end

        %%%%%% Functions used by the Dummy Sampler and both of the SMC
        %%%%%% samplers
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

        function theta_ = SamplePrior(obj, N_samps)
            % samples parameters from the provided priors
            theta_ = zeros(N_samps,obj.P_parameters);
            if isempty(obj.maps)
                for i=1:obj.P_parameters
                    theta_(:,i) = obj.priors{i}.random([N_samps,1]);
                end
            else
                for i=1:size(obj.maps,1)
                    prior_ind = obj.maps{i,1};
                    prior_samps = obj.priors{prior_ind}.random([N_samps,1]);
                    particle_inds = obj.maps{i,2};
                    for j=1:length(particle_inds)
                        theta_(:,particle_inds(j)) = prior_samps(:,j);
                    end
                end
            end

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
                % if this is the finalized distribution, the data is
                % sampled from the kernel density estimation of the
                % posterior (as is used in the pdf calculation)
                samps = cur_particles + mvnrnd(0*cur_particles,diag(obj.bw.^2),len);
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

        function p_ = marginal_pdf(obj, particles, index)
            %{
                Estimates the marginal pdf of parameter at position index 
                using the posterior particles and kernel density
                estimation.
            %}
            p_ = ksdensity(obj.posterior(:,index),particles,'Bandwidth',obj.bw(index));
            
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
            fprintf("~~~~ Posterior Summary (median [95%% CI]) ~~~~\n")
            for p=1:obj.P_parameters
                summary(p,1) = median(obj.posterior(:,p));
                summary(p,2:3) = quantile(obj.posterior(:,p),[0.025,0.975]);
                fprintf("\t%s: %.3f [ %.3f,%.3f ] \n",obj.parameter_labels{p},summary(p,1),summary(p,2),summary(p,3))
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

        function fig = PlotKDEstimate(obj)
            
            if ~obj.calibrated
                error("Posterior must be calibrated before you can plot the convergence")
            end
            
            % for each of the parameters, plot the histograms of particles
            figs = cell(obj.P_parameters,1);
            for i=1:obj.P_parameters
                fi = figure("Position",[100,100,600,300],'Color','w'); hold all
                range = [min(obj.posterior(:,i)), max(obj.posterior(:,i))];
                x_ = linspace(range(1) - 0.1*diff(range), range(2) + 0.1*diff(range), 300);
                y_ = obj.marginal_pdf(x_, i);

                F_ = cumtrapz(x_,y_); % approximate CDF
                iF_ = @(p) fminsearch(@(x) (p - interp1(x_(1:end), F_, x)).^2, mean(x_));
                p_0_025 = iF_(0.025);
                p_0_50 = iF_(0.5);
                p_0_975 = iF_(0.975);

                fill([x_,x_(1)],[y_,y_(1)],'b','FaceAlpha',0.4,'EdgeColor','none');
                errorbar(p_0_50,1.1*max(y_), p_0_025-p_0_50,p_0_975-p_0_50,'horizontal','sk','MarkerFaceColor','k','capsize',0,'linewidth',1)
                plot(x_,y_,'-','Color',[0.2,0.4,0.9],'LineWidth',1.5);
                xlabel(obj.parameter_labels(i));
                ylabel("\pi("+obj.parameter_labels(i)+")")
                set(gca,"FontName","calibri","FontSize",12);
                ylim([0,1.15*max(y_)])
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