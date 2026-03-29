# Using Disasters to Estimate the Impact of Uncertainty

Python replication of Baker, Bloom, and Terry (2024), *Review of Economic Studies* 91, 720–747.

## Overview

This project provides a complete Python reproduction of the original codebase (Stata + MATLAB + Fortran), replicating all tables and figures from the paper bit-for-bit.

### Modules

| Module | Original | Python | Output |
|--------|----------|--------|--------|
| **IV** | Stata | `src/iv/` | Tables 1–6 |
| **IV_VAR** | MATLAB | `src/iv_var/` | Figures 6–7 |
| **LMN_VAR** | Stata + MATLAB | `src/lmn_var/` | Figures 3–5 |
| **MODEL** | Fortran 90 | `src/model/` | Figure 8 |

## Setup

```bash
# Create virtual environment
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

### IV Panel Regressions (Tables 1–6)
```bash
python -m src.iv.panel_iv
```

### IV-VAR Estimation (Figures 6–7)
```bash
python -m src.iv_var.run_estimation
```

### LMN VAR (Figures 3–5)
```bash
python -m src.lmn_var.run_estimation
```

### Model Simulation (Figure 8)
```bash
python -m src.model.run_simulation
```

## Project Structure

```
├── data/
│   ├── IV/              # Panel IV data (.dta)
│   ├── IV_VAR/          # IV-VAR data (.csv)
│   ├── LMN_VAR/         # LMN VAR data (.dta)
│   └── MODEL/
├── src/
│   ├── iv/              # Panel IV regressions
│   ├── iv_var/          # IV-VAR estimation
│   ├── lmn_var/         # Disaster event restrictions VAR
│   ├── model/           # Micro-macro model
│   └── utils/           # Shared utilities
├── original_code/       # Original code (Stata/MATLAB/Fortran)
│   ├── IV/              # Original Stata code for Panel IV
│   ├── IV_VAR/          # Original MATLAB code for IV-VAR
│   ├── LMN_VAR/         # Original Stata + MATLAB code for LMN VAR
│   ├── MODEL/           # Original Fortran 90 code for model
│   └── README.pdf       # Original documentation
├── output/
│   ├── tables/
│   └── figures/
├── tests/
├── requirements.txt
└── README.md
```

## Original Code

The `original_code/` directory contains the original replication code from Baker, Bloom, and Terry (2024). This includes:

| Directory | Language | Description |
|-----------|----------|-------------|
| `IV/` | Stata (`.do`) | Panel IV regression code |
| `IV_VAR/` | MATLAB (`.m`) | Instrumental Variable VAR estimation |
| `LMN_VAR/` | Stata + MATLAB | Local Projections VAR with disaster restrictions |
| `MODEL/` | Fortran 90 (`.f90`) | Micro-macro uncertainty model |

### Original Code Structure

```
original_code/
├── README.pdf                    # Original documentation
├── IV/
│   ├── Panel IV Code.do         # Stata code for panel IV
│   ├── panel_iv_data.dta        # Data file
│   └── dstats.csv               # Summary statistics
├── IV_VAR/
│   ├── STEP1_ESTIMATION.m       # Main estimation script
│   ├── STEP2_GRAPHS.m           # Figure generation
│   ├── VARdata.csv              # Data file
│   ├── BASELINE/                # Baseline specification
│   ├── BOOT/                    # Bootstrap version
│   ├── EARLY/                   # Early sample
│   ├── LATE/                    # Late sample
│   ├── FEWER_LAGS/              # Sensitivity: fewer lags
│   ├── MORE_LAGS/               # Sensitivity: more lags
│   ├── NO_COUNTRY_FE/           # No country fixed effects
│   └── NO_TIME_FE/              # No time fixed effects
├── LMN_VAR/
│   ├── STEP1_STATA_ESTIMATION.do    # Stata OLS estimation
│   ├── STEP2_MATLAB_ESTIMATION.m    # MATLAB VAR estimation
│   ├── STEP3_GRAPHS.m               # Figure generation
│   ├── Dates_and_Data.dta           # Data file
│   ├── BASELINE/                    # Baseline specification
│   ├── 10LAGS/                      # 10 lags version
│   ├── 14LAGS/                      # 14 lags version
│   ├── LOOSER/                      # Looser restrictions
│   ├── TIGHTER/                     # Tighter restrictions
│   ├── NO_COUNTRY_FE/               # No country fixed effects
│   ├── NO_TIME_FE/                  # No time fixed effects
│   ├── REV_ONLY/                    # Revolutions only
│   └── REV_COUP_ONLY/               # Revolutions & coups only
└── MODEL/
    ├── FIRST_STAGE.m            # First stage estimation
    ├── VOL_GROWTH_wrapper.f90   # Fortran model wrapper
    ├── base_lib.f90             # Base library
    └── compile_script.sh        # Compilation script
```

## Reference

Baker, S. R., Bloom, N., & Terry, S. J. (2024). Using Disasters to Estimate the Impact of Uncertainty. *The Review of Economic Studies*, 91(2), 720–747. https://doi.org/10.1093/restud/rdad036

## License

This replication code is provided for academic research purposes. Original code © Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
