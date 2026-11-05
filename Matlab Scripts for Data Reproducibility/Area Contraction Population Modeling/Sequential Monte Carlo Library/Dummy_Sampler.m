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

        transforms      % transformer classes that allows for bounded variables to be transformed to the full real number line.
                        % It also calculates the necessary corrections to the prior probability in the MH step

        parameter_labels    % array of string labels for each the parameters. Used for labeling plots

        priors          % set of priors for the parameters, must be a list of distributions of the same length as your list of scalar vectors
        report_prior    % if a set of priors is provided, then instead of estimating the pdf based on the particles, the prior probability is reported
        
        bw  % bw for the kernel density estimation of the posterior based on Silverman's Rule of Thumb
        calibrated      % boolean variable set to true once GeneratePosterior has been run

    end

    methods
        function obj = Dummy_Sampler(posterior,parameter_labels,maps,transforms,priors,report_prior)
            arguments 
                posterior 
                parameter_labels
                maps = {};
                transforms = {};
                priors = {};
                report_prior = 0;
            end

            obj.posterior = posterior;

            obj.parameter_labels = parameter_labels;

            obj.N_particles = size(posterior,1);
            obj.P_parameters = size(posterior,2);

            obj.maps = maps;
            obj.transforms = transforms;
            obj.priors = priors;
            obj.report_prior = report_prior & ~isempty(priors);

            obj.bw = obj.CalculateBandwidth();
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
                perturbation = mvnrnd(0*cur_particles,diag(obj.bw.^2),len);
                samps = cur_particles + perturbation;
                if ~isempty(obj.transforms)
                    for j=1:length(obj.transforms)
                        ind = obj.transforms{j}{1};
                        transformed_proposal = obj.transforms{j}{2}.transform(cur_particles(:,ind)) + perturbation(:,ind); % transform the particle and add the transformed perturbation
                        samps(:,ind) = obj.transforms{j}{2}.inv_transform(transformed_proposal);
                    end
                end
            end
        end

        function bw = CalculateBandwidth(obj,particles,indices)
            arguments
                obj 
                particles = obj.posterior;
                indices = 1:size(particles,2);
            end

            bw = zeros(length(indices),1);

            bw_scale = (4 / ((obj.P_parameters + 2)*obj.N_particles) )^( 1/(obj.P_parameters + 4) );
            if ~isempty(obj.transforms)
                for j=1:length(obj.transforms)
                    ind = obj.transforms{j}{1};
                    ind_in_vec = find(indices==ind);
                    if ~isempty(ind_in_vec)
                        particles(:,ind_in_vec) = obj.transforms{j}{2}.transform(particles(:,ind_in_vec));
                    end
                end
            end                
            for i=1:length(indices)
                bw(i) = std(particles(:,i))*bw_scale;
            end
        end

        function p_ = pdf(obj,particles,sample_points)
            arguments
                obj 
                particles 
                sample_points = obj.posterior;
            end
            %{
                Estimate the pdf using posterior particles and kernel
                density estimation. Uses priors to ensure that disallowed
                particles recieve 0 pdf.
            %}
            
            bw_i = obj.CalculateBandwidth(sample_points);

            p_prior = obj.CalculatePrior(particles);
            valid_points = p_prior>0;
            if obj.report_prior
                % report the total prior probability for multistep
                % estimation
                p_ = p_prior;
            else
                % estimate the posterior pdf using kernel density
                % estimation
                if isempty(obj.transforms)
                    p_ = mvksdensity(sample_points,particles(valid_points,:),'Bandwidth',bw_i);
                else
                    transformed_sample_points = sample_points;
                    transformed_particles = particles;
                    log_jacobian = 0;

                    for j=1:length(obj.transforms)
                        ind = obj.transforms{j}{1};
                        transformed_sample_points(:,ind) = obj.transforms{j}{2}.transform(sample_points(:,ind));
                        transformed_particles(valid_points,ind) = obj.transforms{j}{2}.transform(particles(:,ind));
                        log_jacobian = log_jacobian + log(obj.transforms{j}{2}.jacobian(particles(valid_points,ind)));
                    end
                    p_transformed = 0*p_prior;
                    p_transformed(valid_points) = log(mvksdensity(transformed_sample_points,transformed_particles(valid_points,:),'Bandwidth',bw_i));                    
                    p_ = exp(p_transformed + log_jacobian);
                end
                    
                is_range = p_prior > 0;
                p_ = p_ .* is_range;
            end
        end

        function p_ = marginal_pdf(obj, indices, particles, sample_points)
            arguments
                obj 
                indices 
                particles 
                sample_points = obj.posterior;
            end
            %{
                Estimates the marginal pdf of parameter at position index 
                (or of parameters at positions indices) using the sample 
                points and kernel density estimation.
            %}
            if size(sample_points,2)>length(indices)
                sample_points = sample_points(:,indices); % pull out the relevant columns
            end

            bw_i = obj.CalculateBandwidth(sample_points,indices);

            log_jacobian = 0;
            if ~isempty(obj.transforms)
                for j=1:length(obj.transforms)
                    ind = obj.transforms{j}{1};
                    ind_in_vec = find(indices==ind);
                    if ~isempty(ind_in_vec)
                        log_jacobian = log_jacobian + log(obj.transforms{j}{2}.jacobian(particles(:,ind_in_vec)));
                        particles(:,ind_in_vec) = obj.transforms{j}{2}.transform(particles(:,ind_in_vec));
                        sample_points(:,ind_in_vec) = obj.transforms{j}{2}.transform(sample_points(:,ind_in_vec));
                    end
                end
            end
            p_ = log(ksdensity(sample_points,particles,'Bandwidth',bw_i)) + log_jacobian;
            p_ = exp(p_);
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

        function fig = PlotKDEstimate(obj,fig,index,sample_points,color,show_95CI)
            arguments
                obj 
                fig 
                index  
                sample_points = obj.posterior;
                color = [0.2,0.4,0.9];
                show_95CI = 1;
            end
            
            if fig.Type=="figure"
                figure(fig); 
            else
                axes(fig);
            end
            hold all

            % pulling out the relevant vector of particles
            if size(sample_points,2)>1
                sample_points = sample_points(:,index);
            end

            % finding the range of points over which to plot
            range = [min(sample_points), max(sample_points)];
            
            % calculating the marginal pdf
            x_ = linspace(range(1) - 0.1*diff(range), range(2) + 0.1*diff(range), 300)';
            % check if any of the extended range is outside of the allowed
            % range given a transformation
            if ~isempty(obj.transforms)
                valid_points = 1:length(x_);
                for j=1:length(obj.transforms)
                    ind = obj.transforms{j}{1};
                    if index==ind
                        x_trans = obj.transforms{j}{2}.transform(x_);
                        valid_points = boolean((imag(x_trans)==0) .* (abs(x_trans) < Inf));
                    end
                end
                x_ = linspace(min(x_(valid_points)), max(x_(valid_points)),300)';
            end
            
            y_ = obj.marginal_pdf(index,x_,sample_points);

            % ensuring we hit zero
            x_ = [x_(1)-mean(diff(x_)); x_; x_(end)+mean(diff(x_))];
            y_ = [0; y_; 0];
            fill([x_;x_(1)],[y_;y_(1)],color,'FaceAlpha',0.8,'EdgeColor','none');
            plot(x_,y_,'-','Color',0.8*color,'LineWidth',1.5);

            if show_95CI
                % calculating the 95% credible interval around the median
                F_ = cumtrapz(x_,y_); % approximate CDF
                iF_ = @(p) fminsearch(@(x) (p - interp1(x_(1:end), F_, x)).^2, mean(x_));
                p_0_025 = iF_(0.025);
                p_0_50 = iF_(0.5);
                p_0_975 = iF_(0.975);
                errorbar(p_0_50,1.1*max(y_), p_0_025-p_0_50,p_0_975-p_0_50,'horizontal','sk','MarkerFaceColor','k','capsize',0,'linewidth',1)
            end
            
            % beautifying the plot
            xlabel(obj.parameter_labels(index));
            ylabel("\pi("+obj.parameter_labels(index)+")")
            set(gca,"FontName","calibri","FontSize",14);
            ylim([0,1.15*max(y_)])
    

        end

        function fig = PlotKDEstimate_All(obj)
            
            if ~obj.calibrated
                error("Posterior must be calibrated before you can plot the convergence")
            end
            
            % for each of the parameters, plot the histograms of particles
            figs = cell(obj.P_parameters,1);
            for i=1:obj.P_parameters
                fi = figure("Position",[100,100,600,300],'Color','w'); hold all
                fi = obj.PlotKDEstimate(fi,i);
                figs{i} = fi;
            end
        end

        function fig = PlotPrior(obj,fig,index,color)
            arguments
                obj 
                fig 
                index 
                color = [0.8,0.1,0.3]; 
            end
            
            prior_ = obj.priors{index};

            x_ = linspace(prior_.icdf(0.01), prior_.icdf(0.99),300)';
            y_ = prior_.pdf(x_);

            if fig.Type=="figure"
                figure(fig); 
            else
                axes(fig);
            end
            hold all
            
            fill([x_(1);x_;x_(end);x_(1)],[0;y_;0;0],color,'FaceAlpha',0.8,'EdgeColor','none');
            plot(x_,y_,'-','Color',0.8*color,'LineWidth',1.5);

            xlabel(obj.parameter_labels(index));
            ylabel("\pi("+obj.parameter_labels(index)+")")
            set(gca,"FontName","calibri","FontSize",14);
            ylim([0,1.15*max(y_)])
            xlim([min(x_),max(x_)])

        end

        function fig = PlotJointPosterior(obj,fig,indices,show_prior)
            arguments
                obj 
                fig 
                indices 
                show_prior = false;
            end
            figure(fig);
            tiledlayout(length(indices),length(indices),"TileSpacing","tight","Padding","compact"); 
            
            for i=1:length(indices)
                ind = sub2ind([length(indices),length(indices)],i,i);
                axi = nexttile(ind);
                % histogram(axi,obj.posterior(:,i),"Normalization","pdf");
                if show_prior
                    obj.PlotPrior(axi,indices(i))
                end
                obj.PlotKDEstimate(axi,indices(i))
                xlim([min(obj.posterior(:,indices(i))),max(obj.posterior(:,indices(i)))])
                ylabel("")
                xlabel("")
                set(gca,'FontSize',12)

                if i==1
                    ylabel(obj.parameter_labels{indices(i)})
                elseif i==length(indices)
                    xlabel(obj.parameter_labels{indices(i)})
                end
                
                for j=i+1:length(indices)
                    ind = sub2ind([length(indices),length(indices)],i,j);
                    axi = nexttile(ind);
                    
                    if show_prior
                        range_x = linspace(obj.priors{indices(i)}.icdf(0.01),obj.priors{indices(i)}.icdf(0.99),300);
                        range_y = linspace(obj.priors{indices(j)}.icdf(0.01),obj.priors{indices(j)}.icdf(0.99),300);
                    else
                        range_x = linspace(min(obj.posterior(:,indices(i))),max(obj.posterior(:,indices(i))),300);
                        range_y = linspace(min(obj.posterior(:,indices(j))),max(obj.posterior(:,indices(j))),300);
                    end
                    [X,Y] = meshgrid(range_x,range_y);
                    F_ = obj.marginal_pdf([indices(i),indices(j)],[X(:),Y(:)]);
                    F_ = reshape(F_,size(X));
                    contourf(X,Y,F_,8,'EdgeColor','none')
                    xlim([range_x(1),range_x(end)])
                    ylim([range_y(1),range_y(end)])
                    if i==1
                        ylabel(obj.parameter_labels{indices(j)})
                    else
                        yticks([])
                    end
                    if j==length(indices)
                        xlabel(obj.parameter_labels{indices(i)})
                    else
                        xticks([])
                    end
                    set(gca,'FontSize',12)
                    colormap('sky')
                    drawnow()
                end
                
            end
        end
    end
end