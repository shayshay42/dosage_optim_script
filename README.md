# Therapeutic Dose Optimization

This repository contains the implementation of a differential equation-based model to optimize therapeutic dosages for a two-drug treatment of glioblastoma.

## Repository Structure
.
├── assets  
│ ├── dosing.jl  
│ ├── model.jl  
│ ├── params.jl  
│ └── utils.jl  
└── src  
  └── optimize.jl  
- Project Root
  - assets
    - dosing.jl
    - model.jl
    - params.jl
    - utils.jl
  - src
    - optimize.jl


## Files and Folders Description

- **assets**: This folder contains the utility scripts necessary for the model. It includes:
  - `dosing.jl`: Handles dosing schedules and dosing amounts for the drug treatment. It uses a spaced list for dosing times and defines dosing functions.
  - `model.jl`: Contains the pharmacokinetic and pharmacodynamic model equations used in the simulation.
  - `params.jl`: Provides the initial parameters for the differential equation model.
  - `utils.jl`: Contains any additional utility functions used in the scripts.

- **src**: This folder includes the primary script that uses the assets to perform the optimization of dosages:
  - `optimize.jl`: Utilizes the utilities from the assets folder and combines them with a solver to find optimal dosage amounts. It defines an objective function, performs sensitivity analysis, and applies an optimization routine.

## How to Run

Run optimize.jl
Caution! the files here are one splitup script and some variables are hard coded in separate files.
