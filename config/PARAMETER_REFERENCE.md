# Configuration Parameter Reference (Paper + Code Mapping)

This file explains **every parameter** in `config/experiment_config.json`:
- what it controls in code
- where it maps to in the paper (if applicable)
- paper/baseline value

Paper used:  
`Using Disasters to Estimate the Impact of Uncertainty` (Review of Economic Studies, 2024)  
local file: `/Users/if/Downloads/Using Disasters to Estimate the Impact of Uncertainty/Using Disasters to Estimate the Impact of Uncertainty.pdf`

## How to run with config

```bash
python run_all.py
python run_all.py --jobs baseline
python run_all.py iv iv_var
python run_all.py --sequential
```

Runner behavior:
- reads `config/experiment_config.json`
- applies `defaults`
- applies per-job overrides
- writes outputs to `output/runs/<job_name>/...`

## Top-level parameters

| Parameter | Meaning | Paper mapping | Paper value | Repo baseline |
|---|---|---|---|---|
| `parallel.enabled` | Enable process-level parallel execution across jobs | Not in paper (engineering) | N/A | `true` |
| `parallel.max_workers` | Max concurrent jobs | Not in paper | N/A | `2` |
| `default_modules` | Modules run when a job does not specify its own module list | Not in paper | N/A | `["iv","iv_var","lmn_var","model"]` |
| `jobs[].name` | Job identifier, also output subfolder name | Not in paper | N/A | user-defined |
| `jobs[].enabled` | Whether job runs by default | Not in paper | N/A | user-defined |

## IV module (`defaults.iv`)

Code path: `src/iv/panel_iv.py`

| Parameter | Meaning | Paper mapping | Paper value | Repo baseline |
|---|---|---|---|---|
| `data_path` | Input Stata panel dataset | Empirical sample for Tables 1-6 | Dataset defined by paper/data packet | `data/IV/panel_iv_data.dta` |
| `standardize_residualized` | Scaling mode for first/second moments | Table notes: first- and second-moment series scaled to residualized unit std over regression sample | Qualitative rule in paper notes | `none` (retain raw coefficient scale in current pipeline) |

Paper text evidence:
- PDF p.10/p.13/p.15-p.17 table notes include:  
  "scaled ... to have residualized unit standard deviation over the regression sample."

## IV_VAR module (`defaults.iv_var`)

Code path: `src/iv_var/estimation.py`

| Parameter | Meaning | Paper mapping | Paper value | Repo baseline |
|---|---|---|---|---|
| `data_path` | Input CSV for IV-VAR | Figure 6-7 sample | Dataset defined by paper/data packet | `data/IV_VAR/VARdata.csv` |
| `n_starts` | Number of deterministic multi-start optimization runs | Not explicit in paper; solver-stability control | N/A | `120` |
| `start_jitter` | Initial-guess perturbation scale for extra starts | Not in paper | N/A | `0.2` |
| `selection_mode` | Candidate selection strategy across multi-start solutions | Paper gives target impact magnitude but not optimizer tie-break | N/A | `paper_anchor` |
| `objective_tie_tol` | Near-optimal objective tolerance for `paper_anchor` | Not in paper | N/A | `2e-12` |
| `target_impact_t1` | Target immediate impact used by `paper_anchor` | Figure 6 text: immediate drop "just over 3.5 percentage points" | ~`-3.5` | `-3.5` |
| `diag_floor` | Minimum diagonal regularity for B22/B33 in candidate filter | Not in paper | N/A | `0.16` |
| `bootstrap_n` | Number of bootstrap replications | In baseline MATLAB code for Figure 6 bands | `150` (code packet baseline) | `150` |
| `bootstrap_block_size` | Mean block length for stationary block bootstrap | In baseline MATLAB code | `25` (code packet baseline) | `25` |

Paper text evidence:
- PDF p.22 (Figure 6 discussion):  
  "immediate drop of just over 3.5 percentage points in GDP growth."

## LMN_VAR module (`defaults.lmn_var`)

Code path: `src/lmn_var/estimation.py`

| Parameter | Meaning | Paper mapping | Paper value | Repo baseline |
|---|---|---|---|---|
| `data_path` | Input Stata dataset for event-restriction VAR | Figures 3-5 | Dataset defined by paper/data packet | `data/LMN_VAR/Dates_and_Data.dta` |
| `n_draws` | Random draws for admissible-set search | Paper discusses admissible-response envelopes and max-G; draw count itself is from code packet | Not explicitly stated in paper text | `1500000` |

Paper text evidence:
- PDF p.20 (Figure 3 discussion): uncertainty shock impact range `1` to `2.5` pp, median around `2%`, and `max G` concept.

Code-packet baseline evidence:
- `original codes and data/BASELINE_var_LMN.m`: `N = 1500000`.

## MODEL module (`defaults.model`)

Code path: `src/model/solve.py`

| Parameter | Meaning | Paper mapping | Paper value | Repo baseline |
|---|---|---|---|---|
| `simplified` | Reduced-grid fast mode | Not a paper parameter; runtime approximation switch | N/A | `false` |
| `do_estimation` | Run structural GMM step | Supplemental MODEL workflow | N/A | `false` |
| `irf_t` | IRF horizon | Figure 8-style plotting horizon | Not explicitly fixed in paper text | `40` |
| `irf_n_sims` | Number of simulated IRF paths to average | Monte Carlo precision/runtime tradeoff | Not explicitly fixed in paper text | `100` |

Important note:
- Code packet README states MODEL is supplemental ("outside published paper" main quantitative tables/figures).

## Job overrides

Any parameter under `jobs[i].<module_name>` overrides `defaults.<module_name>` for that job only.

Example:
- `jobs[i].iv_var.bootstrap_n = 30` runs a faster bootstrap for that job.
- `jobs[i].lmn_var.n_draws = 5000` runs faster admissible-set checks.

## Recommended baseline for paper-facing replication

- `iv_var.target_impact_t1 = -3.5`
- `iv_var.selection_mode = "paper_anchor"`
- `iv_var.bootstrap_n = 150`
- `lmn_var.n_draws = 1500000`
- `model.simplified = false`

