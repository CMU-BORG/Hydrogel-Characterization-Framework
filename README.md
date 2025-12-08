# Hydrogel-Characterization-Framework
Experimental and computational framework for alginate hydrogel characterization using MTS testing and ABC-SMC modeling

# Alginate Hydrogel Characterization Framework

This repository provides an **experimental and computational framework** for characterizing alginate hydrogels using CAD design, 3D printing, MTS mechanical testing, and Approximate Bayesian Computation – Sequential Monte Carlo (ABC-SMC) modeling. It is designed to be reproducible and extensible, enabling future researchers to replicate experiments and build upon the computational methods.

-----------------------------------------------------------------

## 📂 Repository Structure



------------------------------------------------------------------


---

## 🧪 Experimental Workflow

1. **CAD & Printing**
   - Slice geometry (e.g., `Part3.stl`) using Slic3r profile (30G needle, 5 mL syringe).
   - Print alginate samples and perform post-printing crosslinking.

2. **Sample Metadata**
   - Record diameter, mass, height in `Inventory_Final.xlsx`.
   - Label samples:  
     - `S` = low temp group, `T` = high temp group  
     - `0–6` = crosslink time  
     - `a–g` = replicate ID  

3. **Video Analysis**
   - Record video during crosslinking, upload to sample folder.
   - Run `GetFrames.m` (uses `TraceVideo.m`) → saves frames in `frames/`.

4. **Area Calculation**
   - Run `CircleFitting.m` (uses `natsort.m`, `ErrEllipse.m`) → calculates cross-sectional area.
   - Record diameters in `Inventory_Final.xlsx`, area time series in `Areas_Final.csv`.

5. **Mechanical Testing (MTS)**
   - Material stiffness was analyzed on a **universal testing system (MTS Systems, Criterion Model C41)** with a **biobath accessory** using a **50N load cell**.
   - The biobath was filled to a fixed height with **2.5% (w/v) CaCl₂ solution** for aquatic testing.
   - The crosshead was equipped with a compression plate.  
   - An initial test measured ambient buoyancy forces.  
   - Crosslinked samples were fixed to the stage using **double-sided tape and water-resistant glue**.  
   - The crosshead plate was lowered until contact with the sample at a **pre-load of 0.10 mN**.  
   - Samples were cyclically compressed **1.5 mm at 0.1 mm/s for 30 cycles** until peak forces converged within <5%.  
   - Raw force vs. displacement data were exported to CSV for MATLAB analysis.

6. **Data Processing in MATLAB**
   - Force vs. displacement data were corrected for ambient noise and buoyancy forces.  
   - Loading curves from the **last 5 cycles** were isolated.  
   - Plastic deformation from pre-straining resulted in a **zero-contact region (~1.1 mm)**; curves were truncated to include only contact regions.  
   - Sample dimensions were corrected before stress (σ, N/mm²) and strain (ε, mm/mm) calculations.  
   - The stress region occurred over ~10% strain.  
   - Effective stiffness modulus (E_eff) was approximated using a **secant linear slope (0.07 < ε < 0.10)** of the converged stress-strain curves.  
   - The mean slope was recorded as the effective stiffness modulus.

7. **Swelling Analysis**
   - Crosslinking time and temperature effects were quantified by tracking **cross-sectional area changes**.  
   - Frames from videos were sampled every **15 min** and imported to MATLAB.  
   - Each frame was converted to an **HSV pixel map**; a mask isolated the blue dye range.  
   - A size threshold filter reduced noise.  
   - The map was converted to a binary image for edge detection.  
   - Average radii from the map center to detected edges were used to calculate cross-sectional area.

---

## 🔬 Bayesian Modeling (ABC-SMC)

The computational framework uses **Approximate Bayesian Computation – Sequential Monte Carlo (ABC-SMC)** to infer posterior distributions of material stiffness parameters from experimental data. This approach allows us to estimate not only the mean stiffness trajectory but also the full predictive spread of possible outcomes.

### Overview of ABC-SMC
- **Goal:** Infer posterior distributions for parameters governing stiffness evolution:
  - `E0` (initial stiffness)
  - `dE` (stiffness gain after crosslinking)
  - `τ` (time constant of crosslinking kinetics)
- **Challenge:** The likelihood function for viscoelastic hydrogel stiffness is intractable.
- **Solution:** Use ABC-SMC, which compares simulated summary statistics to experimental ones without requiring an explicit likelihood.

### Workflow
1. **Initialization**
   - Define prior distributions for `E0`, `dE`, and `τ`.
   - Set tolerance thresholds for acceptance.

2. **Simulation**
   - For each particle:
     - Sample candidate parameters from priors or proposal kernels.
     - Generate simulated stiffness curves:
       

\[
       E(t) = E_0 + dE \cdot (1 - e^{-t/\tau})
       \]


     - Compute summary statistics:
       - `sum_met_central`: median stiffness trajectory
       - `sum_met_spread`: interquartile range (IQR) of stiffness across replicates

3. **Comparison**
   - Compare simulated summary statistics (`sum_met_central_sim`, `sum_met_spread_sim`) to experimental ones (`sum_met_central`, `sum_met_spread`).
   - Accept particles if the distance between simulated and experimental metrics is below tolerance.

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

