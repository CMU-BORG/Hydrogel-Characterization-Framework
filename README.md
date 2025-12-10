# Hydrogel-Characterization-Framework
Experimental and computational framework for alginate hydrogel characterization using MTS testing and ABC-SMC modeling.

This repository provides an **experimental and computational framework** for characterizing organic hydrogels with high variability. The workflow relies upon material data from MTS mechanical testing, and utilizes Approximate Bayesian Computation – Sequential Monte Carlo (ABC-SMC) modeling to simulate large populations from small raw datasets. This enables robust statistics to be performed and predictions about time-dependent trends to be analyzed on datasets that otherwise are statistically weak at the given sample size. It is designed to be reproducible and extensible, enabling future researchers to replicate experiments and build upon the computational methods for the characterization of additional hydrogel materials.

-----------------------------------------------------------------

## 📂 Repository Structure



------------------------------------------------------------------


---

## 🧪 Experimental Workflow

1. **CAD & Printing**
   - Slice target part geometry (e.g., `Part3.stl`) using Slic3r configuration profile (30G needle, 5 mL syringe).
   - FRESH-print part out of desired hydrogel material, and perform any post-processing steps to bring it to the "uncrosslinked" reference state. 

2. **Sample Metadata**
   - Measure and record initial sample diameter, mass, height in `Inventory_Final.xlsx`.
   - Label samples and create unique sample folder according to the naming structure:  
     - `S` = low temp group, `T` = high temp group  
     - `0–6` = crosslink time  
     - `a–g` = replicate ID

3. **Crosslinking Video Analysis**
   - Crosslinking should be performed according to protocols for the hydrogel material.
   - Record video during crosslinking such that the sample is in focus and appears fully in the frame, save to unique sample folder as MP4.
   - Run `GetFrames.m` (uses `TraceVideo.m`) → saves sample frames from every 15 minutes to `frames/`.

4. **Area Calculation and Swelling Analysis**
   - Run `CircleFitting.m` (uses `natsort.m`, `ErrEllipse.m`) → calculates cross-sectional area in 15 min intervals.
      - Frames from videos are sampled every **15 min** and imported to MATLAB.  
      - Each frame is converted to an **HSV pixel map**; a mask isolates the blue dye range.  
      - A size threshold filter reduces noise.  
      - The map is converted to a binary image for edge detection.  
      - Average radii from the map center to detected edges are used to calculate cross-sectional area.
   - Record cross-sectional areas in `Inventory_Final.xlsx`, area time series in `Areas_Final.csv`.

5. **Mechanical Testing (MTS)**
   - Collect material mechanics data of a sample population
      - Material stiffness can be analyzed on a **universal testing system (MTS Systems, Criterion Model C41)** (or equivalent system) with a **biobath accessory** using a **50N load cell** (or equivalent magnitude load cell with sufficient precision).
      - The biobath should be filled to a fixed, repeatable height for aquatic testing.
      - The crosshead should be equipped with a compression plate.  
      - Perform an initial test to measure ambient buoyancy forces.  
      - Fix crosslinked samples to the stage using **double-sided tape and water-resistant glue**.  
      - Lower the crosshead plate until contact with the sample at a **pre-load of 0.10 mN**.  
      - Cyclically compress samples **1.5 mm (or chosen deformation range) at 0.1 mm/s (or slower to allow for full relaxation) for 30 cycles** until peak forces converged within <5%.  
      - Raw force vs. displacement data are exported to CSV for MATLAB analysis.

6. **Stiffness Data Processing in MATLAB**
   - Run `StiffnessTesting.m` → calculates the effective elastic modulus of tested hydrogel samples.
      - Force vs. displacement data are corrected for ambient noise and buoyancy forces.  
      - Loading curves from the **last 5 cycles** are isolated.  
      - Plastic deformation from pre-straining can result in a **zero-contact region**; Truncate curves to include only contact regions.  
      - Sample dimensions were corrected before stress (σ, N/mm²) and strain (ε, mm/mm) calculations.  
      - The stress region occurred over ~10% strain.  
      - Effective stiffness modulus (E_eff) was approximated using a **secant linear slope (0.07 < ε < 0.10)** of the converged stress-strain curves.  
      - The mean slope was recorded as the effective stiffness modulus.


---

## 🔬 Bayesian Modeling (ABC-SMC)

The computational framework uses **Approximate Bayesian Computation – Sequential Monte Carlo (ABC-SMC)** to infer posterior distributions of material stiffness parameters from experimental data. This approach allows us to estimate not only the mean stiffness trajectory but also the full predictive spread of possible outcomes.

### Overview of ABC-SMC
- **Goal:** Infer posterior distributions for parameters governing stiffness evolution:
  - `E0` (initial stiffness)
  - `dE` (stiffness gain after crosslinking)
  - `τ` (time constant of crosslinking kinetics)
- **Challenge:** The likelihood function for viscoelastic hydrogel stiffness is intractable. The sample size is small, and the variability high, so fit confidence and power is low.
- **Solution:** Use ABC-SMC, which compares simulated summary statistics to experimental ones without requiring an explicit likelihood, to simulate a larger population that can be fit and modeled with high power and confidence.

### Workflow
1. **Initialization**
   - Define prior distributions for `E0`, `dE`, and `τ`.

         mu_E0_mean = ; % [kPa]
         mu_E0_std = ;  % [kPa]
         mu_E0_prior = ;
         
         mu_dE_mean = ; % [kPa]
         mu_dE_std = ;  % [kPa]
         mu_dE_prior = ;
         
         sig_E_lower = ;   % [kPa]   
         sig_E_upper = ;   % [kPa]
         sig_E_prior = ;
         
         rho_E_prior = ;   % [  ]
         
         mean_tau_prior = ;   % [minutes]
         sig_tau_prior = ;    % [minutes]
     
   - Set proposal widths for parameter walks.
     
         meanE_proposal_width = 0.3;   % [kPa]
         stdE_proposal_width = 0.4;    % [kPa]
         rho_proposal_width = 0.1;   % [  ]
         meanTau_proposal_width = 20; % [minutes]
         sigTau_proposal_width = 1;  % [minutes] 

   - Set the number of iterations until convergence.
   - Set tolerance thresholds for acceptance.
      - `epsilon_sum_met_central`: median threshold
      - `epsilon_sum_met_spread`: spread threshold
    


2. **Simulation**
   - For each particle:
     - Sample candidate parameters from priors or proposal kernels.
     - Generate simulated stiffness curves:
       

            E(t) = E_0 + dE \cdot (1 - e^{-t/\tau})


   - Compute summary statistics:
      - `sum_met_central`: median stiffness trajectory
      - `sum_met_spread`: interquartile range (IQR) of stiffness across replicates

3. **Comparison**
   - Compare simulated summary statistics (`sum_met_central_sim`, `sum_met_spread_sim`) to experimental ones (`sum_met_central`, `sum_met_spread`).
   - Accept particles if the distance between simulated and experimental metrics is below tolerance (epsilon).

4. **Resampling**
   - Weight accepted particles by their distance.
   - Resample particles to form the next generation.
   - Proposal kernels are tuned iteratively to improve posterior contraction and acceptance rates.

5. **Posterior Updating**
   - After multiple generations, the particle population approximates the posterior distribution of parameters.
   - Posterior predictive simulations are used to generate confidence intervals for stiffness.

### Key Notes
- **Two levels of variability are modeled:**
  - *Across particles:* posterior uncertainty in mean stiffness trajectory.
  - *Within particles:* spread of possible stiffness values given `sig_E0`, `sig_dE`, and correlation `ρ`.
- **IQR scaling:**  
  - The spread metric (`sum_met_spread`) is based on the interquartile range of stiffness.  
  - .......
- **Final visualization:**  
  - Posterior mean stiffness curve is plotted.
  - Raw data is plotted.
  - Shaded band shows the 95% credible interval derived from the final sampled population, interpolated over a dense time grid (`time_vector = linspace(0,360,360)`).

### MATLAB Implementation
The main script is `Hydrogel_ABC_SMC_normal_E0_dE.m`. Alternative versions of this script include: 
.....



---

## ⚙️ Dependencies

- **MATLAB** (tested with R2021b and later)
- Required MATLAB toolboxes:
  - Statistics and Machine Learning Toolbox

