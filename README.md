# Ordinal Diagnostics - Bachelor Thesis Code

Code repository for the bachelor thesis "Investigating Diagnostic Tools for Ordinal Models"

## Repository Structure

- **`simulation_study/`** - R scripts for the simulation study with different assumption violation scenarios
  - `simulate_CLM.R` - **Main file to run** - executes all simulation scenarios
  - Other files - Supporting functions for the simulation study
- **`empirical_analysis/`** - R scripts for empirical analysis on housing satisfaction and red wine quality datasets
  - `empirical_housing.R` - **Main file** - complete self-contained analysis of housing satisfaction data
  - `empirical_wine.R` - **Main file** - complete self-contained analysis of red wine quality data


## Key Datasets for empirical analysis

- Housing satisfaction survey (MASS::housing)
- Red wine quality dataset

## Main Diagnostic Tools Evaluated

- Surrogate residuals (sure package)
- Proportional odds tests
- Influence and leverage diagnostics
- DHARMa residual diagnostics

All results are fully reproducible with the provided code.
