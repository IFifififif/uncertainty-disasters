# Module/Function Compare (Fresh)

## A. Direct module mapping

### src/iv/panel_iv.py <-> original codes and data/panel_iv.py
- similarity: `0.8377`
- lines(new/old): `927/814`
- common symbols: `16`
- new-only symbols (2): `def _compute_residualized_std, def _scale_to_residualized_unit_std`
- old-only symbols (0): ``

### src/iv_var/estimation.py <-> original codes and data/iv_var_estimation.py
- similarity: `0.3924`
- lines(new/old): `990/611`
- common symbols: `14`
- new-only symbols (4): `def _prepare_panel_data, def _run_spec, def _solve_gmm, def create_mt19937_rng`
- old-only symbols (0): ``

### src/lmn_var/estimation.py <-> original codes and data/estimation.py
- similarity: `0.2926`
- lines(new/old): `723/406`
- common symbols: `9`
- new-only symbols (2): `def _build_companion_form, def create_mt19937_rng`
- old-only symbols (0): ``

### src/utils/regression.py <-> original codes and data/regression.py
- similarity: `0.2832`
- lines(new/old): `751/388`
- common symbols: `6`
- new-only symbols (4): `def areg_ols, def compute_kp_rk_wald_f, def ols_with_classical_se, def partial_out_fwl`
- old-only symbols (0): ``

### src/model/solve.py <-> original codes and data/solve.py
- similarity: `0.1391`
- lines(new/old): `332/302`
- common symbols: `4`
- new-only symbols (6): `def build, def compute_irf, def estimate, def quick_test, def simulate, def solve`
- old-only symbols (6): `def _adjustment_cost, def _bellman_rhs, def _build_grids, def _firm_profit, def simulate_irf, def solve_value_function`

## B. Every src file: best original match

| src file | best old file | similarity | src symbols |
|---|---|---:|---|
| src/__init__.py | original codes and data/src___init__.py | 1.0000 | 0 |
| src/iv/__init__.py | original codes and data/iv___init__.py | 1.0000 | 0 |
| src/iv/__main__.py | original codes and data/__main__.py | 1.0000 | 0 |
| src/iv/panel_iv.py | original codes and data/panel_iv.py | 0.8377 | 18 |
| src/iv_var/__init__.py | original codes and data/iv_var___init__.py | 1.0000 | 0 |
| src/iv_var/__main__.py | original codes and data/iv_var___main__.py | 1.0000 | 0 |
| src/iv_var/estimation.py | original codes and data/iv_var_estimation.py | 0.3924 | 18 |
| src/lmn_var/__init__.py | original codes and data/lmn_var___init__.py | 1.0000 | 0 |
| src/lmn_var/__main__.py | original codes and data/lmn_var___main__.py | 0.7590 | 0 |
| src/lmn_var/estimation.py | original codes and data/estimation.py | 0.2926 | 11 |
| src/model/__init__.py | original codes and data/__main__.py | 0.0364 | 0 |
| src/model/__main__.py | original codes and data/model___main__.py | 0.4413 | 0 |
| src/model/adjustment.py | original codes and data/iv_var_estimation.py | 0.0477 | 12 |
| src/model/ge_solver.py | original codes and data/iv_var_estimation.py | 0.0445 | 12 |
| src/model/gmm.py | original codes and data/iv_var_estimation.py | 0.0312 | 8 |
| src/model/grids.py | original codes and data/iv_var_estimation.py | 0.0445 | 6 |
| src/model/irf.py | original codes and data/regression.py | 0.0380 | 7 |
| src/model/iv_regression.py | original codes and data/regression.py | 0.0547 | 8 |
| src/model/optimizer.py | original codes and data/regression.py | 0.0357 | 5 |
| src/model/params.py | original codes and data/solve.py | 0.0300 | 7 |
| src/model/simulation.py | original codes and data/solve.py | 0.0404 | 13 |
| src/model/solve.py | original codes and data/solve.py | 0.1391 | 10 |
| src/model/vfi.py | original codes and data/solve.py | 0.0350 | 9 |
| src/utils/__init__.py | original codes and data/utils___init__.py | 1.0000 | 0 |
| src/utils/regression.py | original codes and data/regression.py | 0.2832 | 10 |

## C. Function inventory by module

### module `iv`
- `src/iv/__main__.py` (0 symbols)
- `src/iv/panel_iv.py` (18 symbols)
  - class PanelIV@L31, def __init__@L49, def _compute_residualized_std@L118, def _get_fe_arrays@L106, def _prepare_iv_regression@L158, def _run_areg@L242, def _run_iv@L212, def _run_iv_weighted@L593, def _save_table_results@L873, def _scale_to_residualized_unit_std@L139, def load_data@L95, def run_all@L911, def table1_dstats@L300, def table2_baseline@L348, def table3_robustness@L470, def table4_weighting@L646, def table5_media_weightings@L694, def table6_alternative_uncertainty@L740

### module `iv_var`
- `src/iv_var/__main__.py` (0 symbols)
- `src/iv_var/estimation.py` (18 symbols)
  - class IVVAR@L44, def __init__@L64, def _build_moment_vector@L222, def _compute_irf@L368, def _gmm_objective@L305, def _initial_params@L412, def _prepare_panel_data@L81, def _run_spec@L924, def _solve_gmm@L438, def _stationary_block_bootstrap@L781, def bootstrap_se@L654, def create_mt19937_rng@L28, def estimate_baseline@L489, def estimate_robustness@L582, def load_data@L206, def plot_figure6@L848, def plot_figure7@L876, def run_all@L946

### module `lmn_var`
- `src/lmn_var/__main__.py` (0 symbols)
- `src/lmn_var/estimation.py` (11 symbols)
  - class LMNVAR@L45, def __init__@L54, def _build_companion_form@L529, def _check_admissibility@L568, def _get_disaster_restrictions@L550, def create_mt19937_rng@L31, def load_data@L66, def run_all@L688, def step1_estimate_var_fe@L74, def step2_admissible_sets@L259, def step3_generate_figures@L600

### module `model`
- `src/model/__main__.py` (0 symbols)
- `src/model/adjustment.py` (12 symbols)
  - class AdjustmentCostCalculator@L335, def __init__@L343, def _build_matrices@L360, def capital_adjustment_cost@L50, def compute_adjustment_costs_grid@L210, def compute_investment_matrix@L266, def compute_output_matrix@L295, def get_capital_adj_cost@L384, def get_labor_adj_cost@L417, def get_period_return@L448, def labor_adjustment_cost@L127, def output@L24
- `src/model/ge_solver.py` (12 symbols)
  - class Aggregates@L41, class DistributionState@L31, class GESolver@L409, def __init__@L420, def _precompute_matrices@L445, def compute_aggregates@L115, def compute_consumption@L222, def compute_excess_demand@L243, def evolve_distribution@L343, def find_market_clearing_price@L264, def initialize_distribution@L56, def run_simulation@L482
- `src/model/gmm.py` (8 symbols)
  - class GMMSolution@L27, def _ols_coef@L234, def _std_pop@L196, def compute_simulated_moments@L38, def estimate_gmm@L460, def gmm_objective@L400, def objective@L518, def simulate_firms_with_disasters@L268
- `src/model/grids.py` (6 symbols)
  - class StateGrids@L25, def build_full_transition_matrix@L327, def build_grids@L56, def build_tauchen_grid_with_sv@L196, def build_tauchen_transition_sv@L265, def get_state_indices@L398
- `src/model/irf.py` (7 symbols)
  - def _capital_ac@L36, def _labor_ac@L44, def _output@L30, def compute_full_irf@L298, def compute_irf_parallel@L144, def compute_model_moments@L397, def simulate_single_firm@L54
- `src/model/iv_regression.py` (8 symbols)
  - class FirstStageResults@L31, class IVResults@L54, class SecondStageResults@L43, def compute_iv_moments@L283, def first_stage_regression@L62, def run_iv_regression@L223, def second_stage_regression@L146, def simulate_iv_data@L329
- `src/model/optimizer.py` (5 symbols)
  - class PSOConfig@L19, class PSOResult@L33, def nelder_mead_2d@L250, def pso_optimize@L44, def pso_optimize_restart@L213
- `src/model/params.py` (7 symbols)
  - class ModelParameters@L14, def __post_init__@L210, def create_params@L279, def get_data_moments@L219, def get_param_bounds@L254, def get_param_vector@L230, def set_param_vector@L243
- `src/model/simulation.py` (13 symbols)
  - class IRFResults@L81, class SimulationResults@L46, def _build_paths@L585, def box_muller_transform@L91, def compute_figure8_irf@L642, def compute_stock_returns@L277, def create_mt19937_rng@L34, def simulate_all_firms@L313, def simulate_firm_exog@L104, def simulate_firms@L656, def simulate_firms_core@L153, def simulate_firms_with_shock@L661, def simulate_irf@L534
- `src/model/solve.py` (10 symbols)
  - class MicroMacroModel@L52, def __init__@L69, def build@L99, def compute_irf@L162, def estimate@L186, def plot_figure8@L210, def quick_test@L322, def run_all@L258, def simulate@L137, def solve@L107
- `src/model/vfi.py` (9 symbols)
  - class ForecastMatrices@L51, class VFISolution@L36, def build_return_matrices@L110, def compute_ev_matrix@L308, def howard_acceleration_step@L211, def initialize_forecast_matrices@L61, def optimization_step@L361, def solve_vfi@L435, def solve_vfi_simplified@L622

### module `utils`
- `src/utils/regression.py` (10 symbols)
  - def areg_ols@L122, def compute_kp_rk_wald_f@L245, def demean_by_group@L21, def demean_multiple_fe@L45, def format_coef_table@L695, def get_cc_yy_cols@L680, def iv2sls_with_cluster_se@L573, def ols_with_classical_se@L530, def ols_with_cluster_se@L450, def partial_out_fwl@L378

