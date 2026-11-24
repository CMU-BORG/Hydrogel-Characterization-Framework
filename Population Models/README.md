# Population Modeling of Hydrogel Properties

This folder contains the class implementation of the Population Model presented in the manuscript, as well as variants needed to implement the Synthetic Likelihood methods. The base class in this folder can be reused to model different material properties using a combined population-dynamical model.

## 📂 Directory Contents:
* [PopulationModel.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Population%20Models/PopulationModel.m): Base class implementing the generalized model. Implements methods for calculating the likelihood (and log likelihood, for use with the SMC sampler), sampling random dynamical parameter sets, and random dynamic trajectories of the parameter of interest. This class is parameter agnostic and can be reused for any parameter of interest following a first order dynamical law. This can also be used as a template for population models of more complicated systems (either in their dynamics or in their parameter distributions) that are still compatible with the [SMC Library](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/tree/main/Sequential%20Monte%20Carlo%20Library)

* [Mass_PopulationModel.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Population%20Models/Mass_PopulationModel.m): Inherits from PopulationModel. Specific implementation of the population model for our change in mass measurements. The log_likelihood function overrides the base class function to implement a synthetic likelihood where the probability of a joint starting mass and mass at a single time $t$ is estimated using a bivariate kernel density estimate.

* [Stiffness_PopulationModel.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Population%20Models/Stiffness_PopulationModel.m): Inherits from PopulationModel. As with Mass_PopulationModel, this class overrides the log_likelihood function with a kernel density estimate based synthetic likelihood.

* [TestingPopulationModel.m](https://github.com/CMU-BORG/Hydrogel-Characterization-Framework/blob/main/Population%20Models/TestingPopulationModel.m): Testing script to validate the functionality of the PopulationModel class. Can also be used to quickly test and visualized population model predictions, given a set of parameters

## 📐 Model Definition
For a measured hydrogel property (e.g., a geometric or mechanical property), we assume that this property changes over time as the material crosslinks as a first order dynamical system. For property $y_i(t)$ measured from sample $i$ at time $t$ and at temperature $T\in[0,1]$ ($0$: 27 $^o$ C, $1$: 37 $^o$ C), the property evolves as
$$\begin{align}
y_i(t|T) = y_{i,0} + \Delta y_i\biggr(1-e^{-t/\tau_i} \biggr)
\end{align}$$
where $y_{i,0}$, $\Delta y_{i}$, and $\tau_i$ are the starting property, the change in property by steady state, and the characteristic time constant of the crosslinking dynamics. To account for sample-to-sample heterogenity (caused by, e.g., variability in biologically derived stock materials), we assume the sample-level dynamical model parameters $\theta_i=\{y_{i,0},\Delta y_i,\tau_i \}$ is drawn form a generative population model. For this model, we assume that $y_{i,0}$ and $\Delta y_i$ are jointly Normally distributed as
$$\begin{align}
\begin{bmatrix} 
    y_{0,i}\\ 
    \Delta y_i 
\end{bmatrix} \Biggr| T_i
    \sim \mathcal{N} \left( 
\begin{bmatrix} 
    \mu_{y_0}\\ 
    \mu_{\Delta y}(T_i) 
\end{bmatrix}, 
\begin{bmatrix} 
    \sigma_{y_0}^{2} & \rho\sigma_{y_0}\sigma_{\Delta y}(T_i)\\ 
    \rho\sigma_{y_0}\sigma_{\Delta y}(T) & \sigma_{\Delta y}(T_i)^2 
\end{bmatrix} \right)
\end{align}$$
and $\tau_i$ is Gamma distributed as
$$\begin{align}
\tau_i\Bigr|T_i\sim \Gamma\left(\mu_{\tau}(T_i),\sigma_{\tau}(T_i)\right)
\end{align}$$
**Note**: Here, the Gamma distribution is parameterized by its mean and standard deviation, rather than typical location and scale parameters. These are related to the mean and standard deviation as $\alpha_\tau = (\mu_\tau/\sigma_\tau)^2$ and $\beta_\tau=\sigma_\tau^2/\mu_\tau$. The model parameters are functions of the crosslinking temperature to account for possible temperature-induced changes in the final steady state parameter and the rate of crosslinking:
$$
\begin{align} 
    \mu_{\Delta y}(T) &= \mu_{\Delta y}\left(1+\alpha_{\mu}T\right) \\ 
    \sigma_{\Delta y}(T) &= \sigma_{\Delta y}\left(1+\alpha_{\sigma}T\right) \\ 
    \mu_{\tau}(T) &= \mu_{\tau}\left(1+\beta_{\mu}T\right) \\ 
    \sigma_{\tau}(T) &= \sigma_{\tau}\left(1+\beta_{\sigma}T\right) 
\end{align} 
$$
where $\alpha_J$ and $\beta_J$ are the fractional changes in the parameter $J$ between the baseline temperature group and the elevated temperature group. This same formulation could be used for continuously varying temperatures as well. For full details, see "Materials and Methods: Population Modeling of Hydrogel Material Properties" in the manuscript.


## ⚙️ Dependencies

- **MATLAB** (tested with R2021b and later)
- Required MATLAB toolboxes:
  - Statistics and Machine Learning Toolbox