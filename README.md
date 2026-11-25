# Hydrogel-Characterization-Framework
Experimental and computational framework for probabilistic characterization of hydrogel materials using population modeling and Sequential Monte Carlo inference.

This repository provides a **generalizable workflow** that combines experimental measurements with probabilistic population modeling. The framework can be applied to both dynamic measurements (e.g., geometric contraction, swelling, degradation, growth, viscoelastic evolution) and discrete endpoint measurements (e.g., elastic modulus, strength, mass change, biochemical assays). By combining first-order constitutive models with Sequential Monte Carlo (SMC) inference, the workflow enables estimation of population-level material behavior, uncertainty, and expected steady-state properties even when datasets are sparse or highly variable.

Although the experimental case study included in this repository focuses on ionically crosslinked alginate hydrogels, the modeling framework was designed to be material-agnostic and can be adapted to other hydrogel systems, soft biomaterials, and biologically derived tissues.

-----------------------------------------------------------------

## 📂 Repository Structure

For reproducibility, download and extract the repository locally before running any MATLAB scripts. Scripts automatically locate required files relative to the repository root.

```bash

├── README.md
│
├── Population Models
│
│   ├──README.md
│   ├──PopulationModel.m
│   ├──Mass_PopulationModel.m
│   ├──Stiffness_PopulationModel.m
│   └──TestingPopulationModel.m
│
├── Sequential Monte Carlo Library
│   │   │
│   │   ├──README.md
│   │   ├──BernoulliModel.m
│   │   ├──Dummy_Sampler.m
│   │   ├──GammaModel.m
│   │   ├──LogitTransform.m
│   │   ├──LogTransform.m
│   │   ├──N_IG_Model.m
│   │   ├──SMC_Sampler.m
│   │   └──TestingSMC.m
│
├── Matlab Scripts for Data Reproducibility
│
│   ├── CircleFitting.m
│   ├── TraceVideo.m
│   ├── StiffnessTesting.m
│   ├── MonomerRatios.m
│   ├── Plotting.m
│   ├── ErrEllipse.m
│   ├── swft.m
│   │
│   ├──Population Modeling Calibration Scripts
│   │   
│   │   ├── Area Results (timestamped)
│   │   │
│   │   ├── Mass Results (timestamped)
│   │   │
│   │   ├── Stiffness Results (timestamped)
│   │   │
│   │   ├──All_Samplers.m
│   │   ├──AreaDynamicsPopulationModel_SMC.m
│   │   ├──AreaModelParameterFitting.m
│   │   ├──PosteriorInference_AreaTimeConstants.m
│   │   ├──MassDynamicsPopulationModel_SL_SMC.m
│   │   ├──StiffnessDynamicsPopulationModel_SL_SMC.m
│   │   ├──ReplottingResults.m
│   │
│   ├── ExperimentalData
│   │
│   │   ├── Inventory Final.xlsx
│   │   ├── Areas Final.csv
│   │   ├── Anorm.csv
│   │   │
│   │   ├── BUOYANCY
│   │   │   └── DAQ.csv
│   │   │
│   │   ├── FIGURES
│   │   ├── HNMR
│   │   │
│   │   └── Sample Folders
│   │       └── Xtn
│   │           ├── DAQ.csv
│   │           ├── frames
│   │           └── figs
│
└── .gitignore

```
------------------------------------------------------------------


---

## 🧪 Experimental Workflow

1. **Fabrication of Hydrogel Samples**
   - Slice target part geometry using Slic3r configuration profile for desired printer.
   - FRESH-print part out of desired hydrogel material, and perform any post-processing steps to bring it to the "uncrosslinkied/baseline crosslinked" reference state. 

2. **Sample Metadata**
   - Measure and record initial sample diameter, mass, height in `Inventory_Final.xlsx`.
   - Label samples and create unique sample folder according to the naming structure:  
     - `S` = low temp group, `T` = high temp group  
     - `0–6` = crosslink time  
     - `a–g` = replicate ID

3. **Crosslinking Video Analysis**
   - Crosslinking should be performed according to protocols for the hydrogel material.
   - Capture periodic images every 15 minutes during crosslinking such that the sample is in focus and appears fully in the frame, saves sample frames from every 15 minutes to `frames/`.

4. **Area Calculation and Swelling Analysis**
   - Run `CircleFitting.m` (uses `natsort.m`, `ErrEllipse.m`) → calculates cross-sectional area in 15 min intervals.
      - Frames from crosslinking are sampled and imported to MATLAB.  
      - Each frame is converted to an **HSV pixel map**; a mask isolates the blue dye range.  
      - A size threshold filter reduces noise.  
      - The map is converted to a binary image for edge detection.  
      - Average radii from the map center to detected edges are used to calculate cross-sectional area.
   - Record cross-sectional areas in `Inventory_Final.xlsx`, area time series in `Areas_Final.csv`, and save figures of visual cross-sectional area for each frame to `figs/`.
   - Optionally generate dynamic cross-sectional area visualization by running `TraceVideo.m.`, which concatenates/transitions through figures from each frame and saves an `MP4`.

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
      - Sample dimensions are corrected before stress (σ, N/mm²) and strain (ε, mm/mm) calculations.  
      - Effective stiffness modulus (E_eff) is approximated using a **secant linear slope** of the converged stress-strain curves.  
      - The mean slope is recorded as the effective stiffness modulus, and written to the inventory table.
    
7. **Plotting and Visualization in MATLAB**
   - Run `PLotting.m` → visualizes the crosslinking dynamics, swelling trends, and effective elastic modulus of tested hydrogel samples across different experimental groups.
      - Outputs include:
      - Effective Elatic Modulus vs Crosslink Time
      - Geometric Swelling vs Crosslink Time
      - Mass Swelling vs Crosslink Time
      - Initial Geometry
      - Normalized Area over Time (groups and means)
        
--------------------------------------------------------------------

## 🖥️Population Modeling Framework (SL-SMC)

The computational framework uses **Synthetic Likelihood – Sequential Monte Carlo (SL-SMC)** to infer posterior distributions of material parameters from experimental data. This approach allows us to model material behavior as a population rather than a single deterministic response.

### Dynamics Datasets
Examples:

  - Cross-sectional area evolution
  - Swelling trajectories
  - Continuous degradation measurements

For these datasets, population parameters are inferred directly using Sequential Monte Carlo (SMC) methods.

### Sparse Endpoint Measurements
Examples:

  - Effective elastic modulus
  - Mass reduction
  - Mechanical testing requiring destructive assays

For these datasets, the framework combines SMC with likelihood-free inference techniques to estimate population-level dynamics from discrete observations.

## Curent Example Applications
The alginate case study included in this repository demonstrates application of the framework to three distinct material metrics:

###Cross-Sectional Area
Used to quantify:

  - Geometric fidelity
  - Contraction dynamics
  - Temperature-dependent dimensional changes

###Effective Elastic Modulus
Used to quantify:

  - Mechanical performance
  - Stiffness development
  - Population-level variability in material strength

### Mass Reduction
Used to quantify:

  - Water expulsion
  - Network contraction
  - Equilibration behavior during crosslinking

Together, these metrics demonstrate how a common population-modeling framework can reveal shared trends across multiple material properties while accounting for biological and manufacturing variability.

## Sequential Monte Carlo Modeling

Detailed descriptions of:

  - Population model definitions
  - Sequential Monte Carlo algorithms
  - Kernel-density likelihood estimation
  - Temperature-effect parameterization
  - Posterior prediction methods

are provided in the dedicated documentation within:

```bash
Modeling Framework/
├── Population Models/
└── Sequential Monte Carlo Library/
```
---------------------------------------------------------

## ⚙️ Dependencies

- **MATLAB** (tested with R2021b and later)
- Required MATLAB toolboxes:
  - Statistics and Machine Learning Toolbox
  - Image Processing Toolbox (recommended for geometric analysis workflows)

