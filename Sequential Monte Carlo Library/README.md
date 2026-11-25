# Sequential Monte Carlo Estimation of Static Model Parameter Posterior

This library implements a Sequential Monte Carlo (SMC) method to estimate the posterior of an abstract statistical model's parameters given a fixed dataset. **NOTE**: this sampler is not currently set up/validated for SMC-based state estimation. For an arbitary statistical model (implemented as a MATLAB class; see below for requirements), a set of data, and a set of priors, this library automates the calibration process and provides tools for visualization, sampling, and inference of the posterior. For an overview of the SMC method used in this library, see [below](#-overview-of-sequential-monte-carlo-methods).

# 📂 Library Contents:
### Files associated with the SMC Samplers
* [Dummy_Sampler.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/Dummy_Sampler.m): Base class for the Sequential Monte Carlo sampler. This class implements generalized functions associated with sampling from and visualizing distributions that are approximated as a set of equally-weighted particles. 
    
    * Inputs:
        
        * **posterior**: $N\times P$ matrix ($N$: number of particles, $P$: dimensionality of the parameters) containing the particles used to estimate the distribution of interest
        * **parameter_labels**: $N\times 1$ array of strings containing the labels for the parameters of each particle. This text can be coded with LaTex, as these strings are used to population plot labels
        * **maps** (optional): $M\times 2$ cell array of parameter maps used connect multiple indices in a parameter to a single joint prior distribution (e.g., a bivariate Normal or a Normal-Inverse-Gamma distribution). Each row of the cell array contains a single map, where the map is a row cell array with the first element giving the index of the **prior** cell array containing the distribution of interest and the second element containing a vector of indices to the particle vector that will receive that outputs of the multidimensional prior.
        * **transforms** (optional): $Q\times 2$ cell array of variable transformations (either LogTransform or LogitTransform at the moment), where each row of the cell contains the index of the particle element being transformed and an instance of either the LogTransform or LogitTransform classes (see below).
        * **priors** (optional): $P\times 1$ cell array of prior distributions (instances of MATLAB distributions or custom implements models, see below)
        * **report_prior** (optional): Boolean variable indicate whether the sampler should return a kernel density estimate of the probability of an input particles (0) or should return the prior probability (1)

    * Functions:

        * **CalculatePrior(particles)**: calculates the prior probability of the supplied particles given the set of prior specified during instantiation
        * **SamplePrior(N_samps)**: returns an $N_{samps}\times P$ matrix of samples from the prior
        * **SystematicResample(weights,M)**: returns $M$ indices (into, e.g., the posterior matrix) that have been chosen using Systematic Resampling with the provided weights
        * **random(size)**: returns size random samples from the posterior
        * **CalculateBandwidth(particles,indices (optional))**: calculates Silverman's Rule of Thumb bandwidth for a kernel density estimate of the particles. indicies allows you to optionally calculate the bandwidth for a subset of the particle vector, whose elements are specified in the vector of indices
        * **pdf(particles,sample_points (optional))**: estimates the probability density of the supplied particles using Kernel Density Estimation. The kernels are centered on and have their bandwidth set by the sample_points. By default, the sample points are the posterior, but a different set of points (or a subset of the posterior) can be supplied instead
        * **marginal_pdf(indices,particles,sample_points (optional))**: returns the marginal pdf of the particle subset indexed by the indicies using the Kernel Density Estimation. Kernels are centered on the sample_points, which default to the class's posterior
        * **PosteriorSummary()**: returns the median and 95\% credible interval for each variable in the posterior. Also returns a list of strings that can be used to generate a text file containing the values
        * **PlotPosteriorHistograms()**: returns a cell array of figures, where each element contains a histogram plot of one element of the posterior. x axis is labeled using the parameter_labels
        * **PlotKDEstimate(fig,index,sample_points (optional), color (optional), show_95CI(optional))**: plots the kernel density estimated marginal pdf of the variable at index in the sample_points (which defaults to the posterior). Color sets the color of the plot, and show_95CI is a Boolean indicating whether or not to show the median and 95\% credible interval as an error bar above the pdf. The plot is added to the figure fig, which then returned
        * **PlotKDEstimate_All()**: returns a cell array of figures, each containing the output of PlotKDEstimate for one element of the posterior
        * **PlotPrior(fig,index,color (optional))**: plots the prior at index in the list of priors in the figure fig. Color sets the color of the plot. returns fig
        * **PlotJointPosterior(fig,indices,show_prior (optional))**: generates a composite grid of figures (tiledlayout) in fig in which the diagonal elements of show the KD estimate pdfs of the posterior particles at indices. The lower off-diagonals show a contour plot of kernel density estimate of the joint marginal pdf of the two variables corresponding to the diagonals that intersect at the grid location (i.e., the joint KDE at location (i,j) show the joint of particle j and particle i).

* [SMC_Sampler.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/SMC_Sampler.m): Inherits from DummySampler. This class implements the Sequential Monte Carlo method to estimate the posterior distribution of a model's parameters, given a fixed dataset and set of parameter priors. It also serves as a sampler for the calibrated posterior, allow for visualization and sampling using the tools of the DummySampler.

    * Inputs: 

        * **N**: number of particles used in the estimation of the posterior. For simple models with low parameter dimensionality, start around 1000; more particles will be needed for larger dimensionality models or for more complicated posteriors
        * **model**: statistical model captured as an instance of a MATLAB class. This class requires the following minimum functions (use exactly these syntax, as the SMC sampler calls them exactly as listed below; see the example statistical model included in the library):
            
            *  ***setParams(obj,theta)***: takes in the paramter vector theta and updates the model's parameter that are stored as class properties. The exact method of storing the model parameters does not matter, as the SMC sampler only need to be able to set them using this function
            * ***log_likelihood(obj,X)***: takes in a data (the form of which will be model dependent and can contain additional information needed to implement synthetic likelihood methods) and returns the model's log likelihood of the dataset. This same function is used for both true likelihoods and synthetic likelihoods. If you are using synthetic likelihoods, implement your simulations and kernel density estimates here. See [Stiffness_PopulationModel.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Population%20Models/Stiffness_PopulationModel.m) for an example of implementing a synthetic likelihood with kernel density estimation
        
        * **priors**: $P\times 1$ cell array of prior distributions (instances of MATLAB distributions or custom implementations models, see below)
        * **propose_mult**: scalar used to scale the emperical covariance during MCMC perturbation steps (see below). This is used to tune the acceptance rate of your perturbation steps to be around 20\% to 40\%
        * **data**: dataset that will be passed to the model's log likelihood function. The structure of this data can be set in a convenient way for your model as the SMC sampler only ever passes the entire variable to the log likelihood function.
        * **gammas**: vector of geometric bridging coefficients (ranging from 0 to 1). This vector sets the tempering schedule of the sample. Start with around 5-10 values in the vector. If the sampler collapses (very low effective sample size), try increasing the number of generations. Alternatively, you can try a nonlinear tempering schedule, where the values for 0 to 1 change according to a nonlinear function
        * **parameter_labels**: $N\times 1$ array of strings containing the labels for the parameters of each particle. This text can be coded with LaTex, as these strings are used to population plot labels
        * **maps** (optional): $M\times 2$ cell array of parameter maps used connect multiple indices in a parameter to a single joint prior distribution (e.g., a bivariate Normal or a Normal-Inverse-Gamma distribution). Each row of the cell array contains a single map, where the map is a row cell array with the first element giving the index of the **prior** cell array containing the distribution of interest and the second element containing a vector of indices to the particle vector that will receive that outputs of the multidimensional prior.
        * **transforms** (optional): $Q\times 2$ cell array of variable transformations (either LogTransform or LogitTransform at the moment), where each row of the cell contains the index of the particle element being transformed and an instance of either the LogTransform or LogitTransform classes (see below).
        * **priors** (optional): $P\times 1$ cell array of prior distributions (instances of MATLAB distributions or custom implements models, see below)
        * **report_prior** (optional): Boolean variable indicate whether the sampler should return a kernel density estimate of the probability of an input particles (0) or should return the prior probability (1)
    
    * Functions:

        * **GeneratePosterior()**: runs the full SMC algorithm and generates the model parameter posterior given the prescribed data (see below). Saves particles and weights at each generates and populated the posterior variable to allow for visualization
        * **PerturbParticles(cur_particles,generation)**: wrapper to run PerturbParticle_i over all cur_particles in generation using MATLAB's parallel computing toolbox
        * **PerturbParticle_i(cur_particle,proposal_widths,model_i,data_,gammas_,prior,transforms)**: perturbs a single particle (cur_particle) using 10 MCMC steps. proposal_widths contains the scaled covariance matrix (or standard deviation) for the perturbation kernel. model_i, data_, gammas_, prior, and transforms are local copies of **model**, **data**, **gammas**, **priors**, and **transforms** that are passed to the function to help minimize the data shuttling overhead in the parallel loop
        * **CalculateWeightUpdates(particles,generation)**: wrapper to run CalculateWeightUpdate_i over all particles in generation using MATLAB's parallel computing toolbox
        * **CalculateWeightUpdate_i(particle_i,model_i,gammas,data,generation)**: calculates the weight update for particle_i at generation. model_i, gammas, data, are local copies of **model**, **gammas**, and **data** pass to minimize overhead
        * **CalculateWeightedCovariance(weights,particles)**: calculates the weighted covariance of the particles using weights
        * **PlotConvergence()**: plots the covergence of each element of the posterior over every generation as a median and 95\% credible interval to show the convergence of the parameter over the sampling routine
    
* [LogTransform.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/LogTransform.m): variable transformer class that can be used with the SMC_Sampler and DummySampler to convert a strictly positive variable into an infinite domain variable. If you are interested in implementing your own transformer class, this class serves as the minimum example of what functions are **required**. Be sure to match names exactly, as these functions are called by the Samplers.

    * Inputs: **None**
    * Functions:

        * **transform(p)**: returns the log transform of the input p
        * **inv_transform(log_p)**: returns the inverse log transform of log_p
        * **jacobian(p)**: returns the Jacobian of the log transform at point p
        * **inv_jacobian(exp_p)**: returns the inverse jacobian at the point exp_p

* [LogitTransform.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/LogitTransform.m): variable transformer class that can be used with the SMC_Sampler and DummySampler to convert a variable bounded on $[a,b]$ into an infinite domain variable. If you are interested in implementing your own transformer class, this class serves as the minimum example of what functions are **required**. Be sure to match names exactly, as these functions are called by the Samplers.

    * Inputs: 

        * **a**: lower bound of the variable
        * **b**: upper bound of the variable
    
    * Function:

        * **transform(p)**: returns the logit transform of the input p
        * **inv_transform(log_p)**: returns the inverse logit transform of log_p
        * **jacobian(p)**: returns the Jacobian of the logit transform at point p
        * **inv_jacobian(exp_p)**: returns the inverse jacobian at the point exp_p 

* [TestingSMC.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/TestingSMC.m): validation script that confirms the functionality of the SMC sampler using single- and multi-variable statistical models. SMC posterior predictions are confirmed against analytical posterior from conjugate priors (see Binomial model example and Gaussian model example) or posteriors that can be reasonably numerically integrated (see Gamma model example). This script also serves as a minimum example of how to use the SMC_Sampler to estimate a parameter posterior.

### Example statistical models for validation and use as priors

The following classes provide examples of how to implement statistical models in the format requried by the SMC_Sampler. These classes can serve either as priors for parameters or can serve as the **model** in the SMC_Sampler. All classes have the following minimum functions (**NOTE**: the names must match exactly to be compatible with the Sampler classes):

* ExampleClass():
    * Functions:
        * **setParams(theta)**: takes in a parameter vector theta and stores the values in class properties (model specific)
        * **log_likelihood(X)**: returns the log likelihood of the dataset X for the particular model. The structure of X can be model specific as the SMC_sampler always provides a full dataset, rather than components. This function can either calculate an analytical log likelihood based on a known closed form pdf or you can return a synthetic likelihood by sampling from your model, calculating observables that align with the dataset X, and then use, e.g., kernel density estimation to estimate the likelihood of the dataset given the simulation results. For synthetic likelihood methods, ensure that sufficient numbers of sythetic data are generated to give a stable estimate of likelihood across multiple calls of the function
        * **rnd(N)**: returns N random samples from your statistical model. Used when this class is used as a prior

The following MATLAB classes all have the above functions.

* [BernoulliModel.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/BernoulliModel.m): class implementation of a Bernoulli random variable model. 
* [GammaModel.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/GammaModel.m): class implementation of a Gamma random variable model. Takes as an optional input a Boolean variable (MS_param; default: False) that, when True, allows you to use a Gamma model that is parameterized by its mean and standard deviation, rather than the standard shape and location parameters.
* [GaussianModel.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/GaussianModel.m): class implementation of a Normal (Gaussian) random variable model. Takes as an optional input a Boolean variable (input_variance) that allows you to parameterize the model by its mean and variance (rather than its standard deviation). This allows the use of a Normal-Inverse-Gamma conjugate prior.
* [N_IG_Model.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Sequential%20Monte%20Carlo%20Library/N_IG_Model.m): class implementation of a Normal-Inverse-Gamma random variable model. This joint distribution serves as a conjugate prior for a Gaussian random variable with unknown mean and variance.


# 📐 Overview of Sequential Monte Carlo Methods
For a complete introduction to the methods used in this library, see the methods and supplemental materials for the paper linked in the root of this repository. For a given statistical model with a probability density function $f(x | \theta)$, parameterized by the parameter vector $\theta$, we wish to estimate the parameter vector $\theta$, and our uncertainty in it, using some collected dataset $X$. Within a Bayesian framework, this comes in the form of estimating the posterior distribution $\pi(\theta | X)$ using Bayes's Theorem:

$$
\begin{align}
\pi(\theta|X) = \frac{f(X|\theta)\pi(\theta)}{\pi(X)}
\end{align}
$$

where $\pi(\theta)$ is the prior distribution on the model parameters, and $\pi(X)$ is the total probability of the data given the model and serves as a normalizing constant. For most statistical models, sampling directly from the posterior is intractable. However, we can sample from the prior or from another distribution and then update the sample until, in aggregate, our samples represent the posterior. In Sequential Monte Carlo methods, we aim to estimate the posterior distribution using a set of weighted particles ($\{w_i,\theta_i\}$). For a properly calibrated set of weights and particles, we can generate new samples from the posterior, estimate probabilities in the posterior, and calculate estimates of functions on the posterior (e.g., the mean or expectation of a function of $\theta$). 

To achieve this, Sequential Monte Carlo breaks the problem into a series of transition steps, where the weights of particles are progressively updated to approximate intermediate distributions ($\pi_g(\theta|X) \propto f(X|\theta)^{\gamma_g}\pi(\theta)$). When the weights of many particles drop too low (causing the effective sample size to drop below a threshold), all particles are resampled according to their weights. This discards low weight particles and replaces them with higher weight ones. Then to prevent many particles from having the same value, a mutation step, consisting of a series of Markov Chain Monte Carlo (MCMC) steps, occurs to remix the particles. 

```
┌──────────────────────────────────────────────────────────────────────────────────────────────┐
│ Algorithm 1: Sequential Monte Carlo Sampling                                                 │
├──────────────────────────────────────────────────────────────────────────────────────────────┤
│ Input: Observed data D, population model M, prior distributions π(θ), proposal widths,       │
│        number of particles N, number of generations G, and a tempering schedule              │
│        {γ_g} ∈ [0,1]                                                                         │
│                                                                                              │
│ Output: Set of equal-weight samples from the posterior π(θ|X)                                │
├──────────────────────────────────────────────────────────────────────────────────────────────┤
│ 1:  Initialize the particle population ({θ_i}_1) by sampling the prior and setting the       │
│     weights to be even (wi = 1/N ∀ i)                                                        │
│                                                                                              │
│ 2:  for i = 2 to G do                                                                        │
│                                                                                              │
│ 3:      Calculate updated weights (w_(i,g) = w_(i,g−1)π_g(θ_i|Y)/π_(g−1)(θ_i|Y)              │
│                                           = w_(g,i−1)f(Y|θ_i)^(γ_g − γ_(g−1)) ∀ i)           │
│                                                                                              │
│ 4:      Calculate effective sample size (ESS = 1/Σ_i w_i²)                                   │
│                                                                                              │ 
│ 5:      if ESS < 0.75N                                                                       │
│                                                                                              │
│ 6:          Systematically resample particles according to w_(i,g)                           │
│                                                                                              │
│ 7:          Reset weights to w_(i,g) = 1/N                                                   │
│                                                                                              │
│ 8:      end if                                                                               │
│                                                                                              │
│ 9:      for j = 1 to 10 do                                                                   │
│                                                                                              │
│ 10:         Perturb particles with mean-zero Gaussian proposals (proposal                    │
│             covariance ∝ weighted covariance of the particles)                               │
│                                                                                              │
│ 11:         Accept particle according to Metropolis-Hastings ratio (see paper)               │
│                                                                                              │
│ 12:     end for                                                                              │
│                                                                                              │
│ 13:     Save intermediate population and weights                                             │
│                                                                                              │
│ 14: end for                                                                                  │ 
│                                                                                              │ 
│ 15: Final systematic resampling to obtain unweighted posterior sample                        │
│                                                                                              │ 
│ 16: return {θ_(i,g)}                                                                         │
└──────────────────────────────────────────────────────────────────────────────────────────────┘
```
