# Using Disasters to Estimate the Impact of Uncertainty

Python replication of Baker, Bloom, and Terry (2024), *Review of Economic Studies* 91, 720–747.

## Overview

This project provides a high-fidelity Python reproduction of the original codebase (Stata + MATLAB + Fortran) and paper results.  
Some modules (especially IV-VAR / LMN-VAR) involve numerically equivalent solution sets, so exact floating-point identity with MATLAB/Stata/Fortran is not guaranteed.

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

### Unified Configuration (recommended)

All runtime parameters are centralized in:

`config/experiment_config.json`

Detailed parameter-by-parameter documentation (paper mapping + baseline values):

`config/PARAMETER_REFERENCE.md`

You can manage in one place:
- data file names/paths for each module (`data_path`)
- estimation/selection parameters (for example IV-VAR multi-start and selection mode)
- LMN draw count
- MODEL runtime options
- multi-dataset jobs and parallel execution (`jobs`, `parallel`)

Run with unified config:

```bash
python run_all.py
```

Optional:

```bash
# Run selected modules for enabled jobs
python run_all.py iv iv_var

# Run only specific jobs
python run_all.py --jobs baseline

# Force sequential run (disable process parallelization for this run)
python run_all.py --sequential
```

### Configuration Notes

- `experiment_config.json` includes inline `__doc__` fields for quick reference.
- For each parameter, `PARAMETER_REFERENCE.md` states:
  - what it controls in code,
  - which paper section/figure/table it corresponds to,
  - paper value (if explicitly provided),
  - current repository baseline value.
- Parameters not explicitly set by the paper are marked as engineering/runtime controls in `PARAMETER_REFERENCE.md`.

### Reproducibility Scope

- `IV` (Tables 1-6): point-estimate replication target.
- `IV_VAR` / `LMN_VAR` (Figures 3-7): set-identification and optimization-path sensitivity exist; use config baselines for stable paper-facing comparisons.
- `MODEL` (Figure 8): supplemental module in the original packet; treat as mechanism/qualitative validation unless explicitly benchmarking a fixed runtime profile.

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
│   ├── figures/
│   └── runs/            # Per-job outputs from unified runner
├── config/
│   └── experiment_config.json
├── tests/
├── requirements.txt
├── run_all.py
└── README.md
```

## Reference

Baker, S. R., Bloom, N., & Terry, S. J. (2024). Using Disasters to Estimate the Impact of Uncertainty. *The Review of Economic Studies*, 91(2), 720–747. https://doi.org/10.1093/restud/rdad036

## License

This replication code is provided for academic research purposes. Original code © Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
