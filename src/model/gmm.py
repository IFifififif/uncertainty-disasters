"""
GMM Estimation for BBT (2024) Model.

This module implements the indirect inference estimation approach:
1. Simulate model with given disaster parameters
2. Compute moments matching the data
3. Compare simulated moments to data moments
4. Optimize parameters via PSO

The model matches 20 moments:
- 4 disaster types x 2 moments (levels, volatility) x 2 samples (macro, micro)
- Plus 2 second-stage coefficients x 2 samples
"""

import numpy as np
from typing import Tuple, Dict, Callable, Optional
from dataclasses import dataclass

from .params import ModelParameters
from .grids import StateGrids, build_grids
from .vfi import VFISolution, solve_vfi_simplified
from .simulation import simulate_firms
from .optimizer import PSOConfig, pso_optimize


@dataclass
class GMMSolution:
    """Container for GMM estimation results."""
    
    x_opt: np.ndarray        # Optimal parameters
    gmm_value: float         # GMM objective at optimum
    converged: bool
    iterations: int
    simulated_moments: np.ndarray
    data_moments: np.ndarray


def compute_simulated_moments(
    params: ModelParameters,
    grids: StateGrids,
    vfi_sol: VFISolution,
    disaster_levels: np.ndarray,
    disaster_unc_probs: np.ndarray,
    T: int = None
) -> np.ndarray:
    """
    Compute simulated moments from model.
    
    Moments organization (20 total):
    1-4:   First stage, levels LHS, macro (4 disaster types)
    5-8:   First stage, vol LHS, macro
    9-10:  Second stage, macro (first moment, second moment)
    11-14: First stage, levels LHS, micro
    15-18: First stage, vol LHS, micro
    19-20: Second stage, micro
    
    Parameters
    ----------
    params : ModelParameters
        Model parameters.
    grids : StateGrids
        State grids.
    vfi_sol : VFISolution
        VFI solution.
    disaster_levels : np.ndarray
        Impact of disasters on levels (4 values).
    disaster_unc_probs : np.ndarray
        Impact of disasters on uncertainty (4 values).
    T : int, optional
        Simulation length.
    
    Returns
    -------
    moments : np.ndarray
        Simulated moments (20 values).
    """
    p = params
    
    if T is None:
        T = p.numdiscard + p.Ncountries * p.Tper
    
    # Run simulation
    sim_results = simulate_firms_with_disasters(
        params, grids, vfi_sol,
        disaster_levels, disaster_unc_probs,
        T=T
    )
    
    # Build moments using the SAME column construction + OLS flow used by
    # Fortran->MATLABdata.csv + FIRST_STAGE.m.
    moments = np.zeros(20)
    eps = 1e-10

    # 1) Build a panel sample analogous to Fortran->MATLABdata bridge.
    n_countries = max(1, int(getattr(params, 'Ncountries', 1)))
    t_per = max(1, int(getattr(params, 'Tper', 1)))
    start = int(getattr(params, 'numdiscard', 0))
    n_keep = n_countries * t_per
    end = min(len(sim_results.Y_sim), start + n_keep)
    if end - start < 32:
        return moments

    y_level = sim_results.Y_sim.astype(np.float64)

    inst_full = sim_results.disaster_indicators
    if inst_full is None:
        instruments_full = np.zeros((len(y_level), 4), dtype=np.float64)
    else:
        instruments_full = inst_full.astype(np.float64)

    # 2) Reconstruct Fortran-style time series used to write MATLABdata.csv.
    n_total = len(y_level)
    growth_sim = np.zeros(n_total)
    first_sim = np.zeros(n_total)
    second_sim = np.zeros(n_total)
    growth_sim_yr = np.zeros(n_total)
    first_sim_yr = np.zeros(n_total)
    second_sim_yr = np.zeros(n_total)

    # Use model-implied public-firm returns when available; fallback to GDP growth proxy.
    if getattr(sim_results, "returnfirm", None) is not None:
        ret = np.asarray(sim_results.returnfirm, dtype=np.float64)
        if ret.ndim == 2 and ret.shape[0] == n_total and ret.shape[1] > 0:
            first_sim = np.nanmean(ret, axis=1)
        else:
            first_sim[1:] = np.log(np.clip(y_level[1:], eps, None) / np.clip(y_level[:-1], eps, None))
    else:
        first_sim[1:] = np.log(np.clip(y_level[1:], eps, None) / np.clip(y_level[:-1], eps, None))

    # Fortran GDP growth: 100*log(Y_t / Y_{t-1})
    growth_sim[1:] = 100.0 * np.log(np.clip(y_level[1:], eps, None) / np.clip(y_level[:-1], eps, None))

    # Fortran second_sim before overwrite (quarterly variance proxy).
    for t in range(3, n_total):
        wnd = first_sim[t - 3:t + 1]
        mu = np.mean(wnd)
        v = np.mean(wnd ** 2) - mu ** 2
        second_sim[t] = np.sqrt(max(v, 0.0))

    # Annualized rolling averages.
    for t in range(4, n_total):
        second_sim_yr[t] = 0.25 * (second_sim[t] + second_sim[t - 1] + second_sim[t - 2] + second_sim[t - 3])
        growth_sim_yr[t] = 0.25 * (growth_sim[t] + growth_sim[t - 1] + growth_sim[t - 2] + growth_sim[t - 3])
        first_sim_yr[t] = 0.25 * (first_sim[t] + first_sim[t - 1] + first_sim[t - 2] + first_sim[t - 3])

        # Fortran overwrite of SECONDsim(t): annual variance around 1% benchmark.
        second_sim[t] = (
            (first_sim[t] - 0.01) ** 2
            + (first_sim[t - 1] - 0.01) ** 2
            + (first_sim[t - 2] - 0.01) ** 2
            + (first_sim[t - 3] - 0.01) ** 2
        )

    # Fortran overwrite of FIRSTsim(t) with annual average.
    first_sim[4:] = first_sim_yr[4:]

    # 3) Build DATAMAT rows exactly like Fortran MATLABdata.csv writer:
    # growthsim at ct, but all RHS and instruments at ct-1.
    rows = []
    for c_idx in range(n_countries):
        for t_idx in range(t_per):
            ct = start + c_idx * t_per + t_idx
            if ct <= 0 or ct >= n_total:
                continue
            lag = ct - 1
            rows.append([
                1.0,
                growth_sim[ct],
                growth_sim_yr[ct],
                first_sim[lag],
                second_sim[lag],
                first_sim_yr[lag],
                second_sim_yr[lag],
                instruments_full[lag, 0] if lag < len(instruments_full) else 0.0,
                instruments_full[lag, 1] if lag < len(instruments_full) else 0.0,
                instruments_full[lag, 2] if lag < len(instruments_full) else 0.0,
                instruments_full[lag, 3] if lag < len(instruments_full) else 0.0,
                float(c_idx + 1),  # country code (1-based)
                float(t_idx + 1),  # within-country time index (1-based)
            ])

    if len(rows) <= 30:
        return moments
    datamat = np.asarray(rows, dtype=np.float64)

    # 4) FIRST_STAGE.m transforms.
    def _std_pop(x: np.ndarray) -> float:
        x = x[np.isfinite(x)]
        if x.size <= 1:
            return np.nan
        return float(np.std(x, ddof=0))

    sd4 = _std_pop(datamat[:, 3])
    if np.isfinite(sd4) and sd4 > eps:
        datamat[:, 3] /= sd4
    sd6 = _std_pop(datamat[:, 5])
    if np.isfinite(sd6) and sd6 > eps:
        datamat[:, 5] /= sd6

    datamat[:, 4] = np.log(np.clip(datamat[:, 4], eps, None))
    sd5 = _std_pop(datamat[:, 4])
    if np.isfinite(sd5) and sd5 > eps:
        datamat[:, 4] /= sd5

    datamat[:, 6] = np.log(np.clip(datamat[:, 6], eps, None))
    sd7 = _std_pop(datamat[:, 6])
    if np.isfinite(sd7) and sd7 > eps:
        datamat[:, 6] /= sd7

    # 5) Dummy matrices like MATLAB dummyvar (full set, no dropped base).
    country = datamat[:, 11].astype(int)
    time = datamat[:, 12].astype(int)
    cc = np.eye(country.max(), dtype=np.float64)[country - 1]
    tt = np.eye(time.max(), dtype=np.float64)[time - 1]

    valid = np.isfinite(datamat[:, [2, 3, 4, 5, 6, 7, 8, 9, 10]]).all(axis=1)
    if valid.sum() <= 20:
        return moments
    dm = datamat[valid]
    ccv = cc[valid]
    ttv = tt[valid]

    Z = np.column_stack([dm[:, 7:11], ccv, ttv])

    def _ols_coef(X: np.ndarray, y: np.ndarray) -> np.ndarray:
        XtX = X.T @ X
        Xty = X.T @ y
        try:
            return np.linalg.solve(XtX, Xty)
        except np.linalg.LinAlgError:
            return np.linalg.pinv(XtX) @ Xty

    # Macro sample
    beta_fm_macro = _ols_coef(Z, dm[:, 3])  # fret
    beta_sm_macro = _ols_coef(Z, dm[:, 4])  # sretcs
    fm_macro_hat = Z @ beta_fm_macro
    sm_macro_hat = Z @ beta_sm_macro
    X2_macro = np.column_stack([fm_macro_hat, sm_macro_hat, ccv, ttv])
    beta2_macro = _ols_coef(X2_macro, dm[:, 2])  # growthsimyr

    # Micro sample
    beta_fm_micro = _ols_coef(Z, dm[:, 5])  # fretann
    beta_sm_micro = _ols_coef(Z, dm[:, 6])  # srettsann
    fm_micro_hat = Z @ beta_fm_micro
    sm_micro_hat = Z @ beta_sm_micro
    X2_micro = np.column_stack([fm_micro_hat, sm_micro_hat, ccv, ttv])
    beta2_micro = _ols_coef(X2_micro, dm[:, 2])  # growthsimyr

    moments[0:4] = beta_fm_macro[:4]
    moments[4:8] = beta_sm_macro[:4]
    moments[8:10] = beta2_macro[:2]
    moments[10:14] = beta_fm_micro[:4]
    moments[14:18] = beta_sm_micro[:4]
    moments[18:20] = beta2_micro[:2]

    return moments


def simulate_firms_with_disasters(
    params: ModelParameters,
    grids: StateGrids,
    vfi_sol: VFISolution,
    disaster_levels: np.ndarray,
    disaster_unc_probs: np.ndarray,
    T: int
):
    """
    Simulate firms with disaster events.
    
    Disaster types:
    0: Natural disaster (probability 0.242)
    1: Political shock (probability 0.03)
    2: Revolution (probability 0.011)
    3: Terrorist attack (probability 0.008)
    
    Each disaster:
    - Reduces GDP level by disaster_levels[i]
    - Increases uncertainty with probability disaster_unc_probs[i]
    """
    p = params
    g = grids
    
    # Disaster probabilities (Fortran DISASTERprobs)
    disaster_probs = np.array([0.242, 0.03, 0.011, 0.008], dtype=np.float64)
    
    # CRITICAL FIX: Use RandomState matching Fortran's random_number
    rng = np.random.RandomState(2501)
    
    # Exogenous states
    a_pos = np.zeros(T, dtype=int)
    s_pos = np.zeros(T, dtype=int)
    disaster_occurred = (rng.random((T, 4)) <= disaster_probs[None, :]).astype(np.float64)
    achg_shocks = rng.random(T)
    schg_shocks = rng.random(T)
    a_shocks = rng.random(T)
    s_shocks = rng.random(T)

    a_pos[0] = max(0, min(p.anum - 1, p.ainit - 1))
    s_pos[0] = max(0, min(p.snum - 1, p.sinit - 1))
    erg_meana = float(np.mean(np.log(np.clip(g.a_grid, 1e-12, None))))
    
    for t in range(1, T):
        # Baseline uncertainty transition (from pr_mat_s), then disaster-induced switch.
        pr_s = g.pr_mat_s[s_pos[t - 1], :]
        s_pos[t] = int(np.searchsorted(np.cumsum(pr_s), s_shocks[t], side="right"))
        s_pos[t] = min(max(s_pos[t], 0), p.snum - 1)

        second_mom_prob = float(np.sum(disaster_occurred[t, :] * disaster_unc_probs))
        second_mom_prob = min(max(second_mom_prob, 0.0), 1.0)
        if schg_shocks[t] <= second_mom_prob:
            s_pos[t] = 1

        # Baseline aggregate productivity transition.
        pr_a = g.pr_mat_a[a_pos[t - 1], :, s_pos[t - 1]]
        a_pos[t] = int(np.searchsorted(np.cumsum(pr_a), a_shocks[t], side="right"))
        a_pos[t] = min(max(a_pos[t], 0), p.anum - 1)

        # Disaster first-moment effect (Fortran uncexogsim logic).
        if disaster_occurred[t, :].sum() > 0.0:
            tot_first_mom = float(p.sigmaa * np.sum(disaster_occurred[t, :] * disaster_levels))
            if tot_first_mom >= 0.0:
                denom = np.log(max(g.a_grid[-1], 1e-12)) - erg_meana
                prob = (tot_first_mom / denom) if abs(denom) > 1e-12 else 0.0
                prob = min(max(prob, 0.0), 1.0)
                if achg_shocks[t] <= prob:
                    a_pos[t] = p.anum - 1
            else:
                denom = np.log(max(g.a_grid[0], 1e-12)) - erg_meana
                prob = (tot_first_mom / denom) if abs(denom) > 1e-12 else 0.0
                prob = min(max(prob, 0.0), 1.0)
                if achg_shocks[t] <= prob:
                    a_pos[t] = 0

    # Run the full simulation stack with exogenous paths overridden by
    # Fortran-style disaster-adjusted (a,s) paths.
    sim_results = simulate_firms(
        params,
        grids,
        vfi_sol,
        T=T,
        seed=2501,
        verbose=False,
        a_pos_override=a_pos,
        s_pos_override=s_pos,
        disaster_indicators=disaster_occurred,
    )
    return sim_results


def gmm_objective(
    x: np.ndarray,
    params: ModelParameters,
    grids: StateGrids,
    data_moments: np.ndarray,
    data_se: np.ndarray,
    vfi_sol: VFISolution = None
) -> float:
    """
    GMM objective function.
    
    GMM = sum_i ((sim_mom_i - data_mom_i) / data_se_i)^2
    
    Parameters
    ----------
    x : np.ndarray
        Parameter vector (8 values):
        - x[0:4]: disaster levels impacts
        - x[4:8]: disaster uncertainty impacts
    params : ModelParameters
        Model parameters.
    grids : StateGrids
        State grids.
    data_moments : np.ndarray
        Data moments to match (20 values).
    data_se : np.ndarray
        Standard errors of data moments.
    vfi_sol : VFISolution, optional
        Pre-computed VFI solution. If None, solves VFI.
    
    Returns
    -------
    gmm_value : float
        GMM objective value (lower is better).
    """
    # Extract parameters
    disaster_levels = x[0:4]
    disaster_unc_probs = x[4:8]
    
    # Check bounds
    if np.any(disaster_unc_probs < 0) or np.any(disaster_unc_probs > 1):
        return 1e10  # Penalty for invalid probabilities
    
    # Solve VFI if not provided
    if vfi_sol is None:
        vfi_sol = solve_vfi_simplified(params, grids, verbose=False)
    
    # Compute simulated moments
    sim_moments = compute_simulated_moments(
        params, grids, vfi_sol,
        disaster_levels, disaster_unc_probs
    )
    
    # GMM objective
    diff = (sim_moments - data_moments) / data_se
    gmm_value = np.sum(diff ** 2)
    
    return gmm_value


def estimate_gmm(
    params: ModelParameters,
    grids: StateGrids = None,
    x_init: np.ndarray = None,
    max_evals: int = 100,
    method: str = "pso",
    pso_max_iter: int = 5000,
    verbose: bool = True
) -> GMMSolution:
    """
    Estimate model parameters via GMM.
    
    Uses PSO by default to match the original Fortran wrapper.
    
    Parameters
    ----------
    params : ModelParameters
        Model parameters.
    grids : StateGrids, optional
        State grids. If None, builds from params.
    x_init : np.ndarray, optional
        Initial parameter guess.
    max_evals : int
        Maximum function evaluations.
    verbose : bool
        Print progress.
    
    Returns
    -------
    GMMSolution
        Estimation results.
    """
    from scipy.optimize import minimize
    
    if grids is None:
        grids = build_grids(params)
    
    # Initial guess
    if x_init is None:
        x_init = params.get_param_vector()
    
    # Data moments and standard errors
    data_moments = params.get_data_moments()
    data_se = params.data_se
    
    # Bounds
    lb, ub = params.get_param_bounds()
    bounds = list(zip(lb, ub))
    
    if verbose:
        print("Starting GMM estimation...")
        print(f"  Initial parameters: {x_init}")
        print(f"  Data moments: {data_moments}")
    
    # Solve VFI once (simplified for speed)
    vfi_sol = solve_vfi_simplified(params, grids, verbose=False)
    
    # Objective function
    def objective(x):
        return gmm_objective(x, params, grids, data_moments, data_se, vfi_sol)
    
    # Optimize
    method_l = method.lower()
    if method_l == "pso":
        # Fortran VOL_GROWTH_wrapper.f90:
        # npart=75, phi=(2.05,2.05), seed=8791, maxpsoit=5000.
        pso_cfg = PSOConfig(
            npart=75,
            max_iter=pso_max_iter,
            x_tol=1e-3,
            f_tol=1e-3,
            x_quick_tol=1e-3,
            x_quick_num=5,
            phi=(2.05, 2.05),
            seed=8791,
        )
        pso_res = pso_optimize(objective, lb, ub, config=pso_cfg, verbose=verbose)
        x_opt = pso_res.x_opt
        f_opt = float(pso_res.f_opt)
        converged = bool(pso_res.converged)
        iterations = int(pso_res.iterations)
    else:
        result = minimize(
            objective,
            x_init,
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': max_evals, 'disp': verbose}
        )
        x_opt = result.x
        f_opt = float(result.fun)
        converged = bool(result.success)
        iterations = int(result.nit)
    
    # Compute final moments
    final_moments = compute_simulated_moments(
        params, grids, vfi_sol,
        x_opt[0:4], x_opt[4:8]
    )
    
    return GMMSolution(
        x_opt=x_opt,
        gmm_value=f_opt,
        converged=converged,
        iterations=iterations,
        simulated_moments=final_moments,
        data_moments=data_moments
    )
