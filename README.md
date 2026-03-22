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
├── output/
│   ├── tables/
│   └── figures/
├── tests/
├── requirements.txt
└── README.md
```

## Reference

Baker, S. R., Bloom, N., & Terry, S. J. (2024). Using Disasters to Estimate the Impact of Uncertainty. *The Review of Economic Studies*, 91(2), 720–747. https://doi.org/10.1093/restud/rdad036

## License

This replication code is provided for academic research purposes. Original code © Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
