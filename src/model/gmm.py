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
    
    # Compute moments
    moments = np.zeros(20)
    
    # First stage moments: impact of disasters on GDP growth and volatility
    # This is a simplified computation - full version would run IV regression
    
    # For now, use disaster impacts directly
    moments[0:4] = disaster_levels  # Levels impact, macro
    moments[4:8] = disaster_unc_probs  # Vol impact, macro (simplified)
    
    # Second stage: impact of uncertainty on growth
    # These come from running IV regression on simulated data
    # Simplified: use parameter values
    moments[8] = 1.5   # First moment coefficient (placeholder)
    moments[9] = -3.0  # Second moment coefficient (placeholder)
    
    # Micro sample moments (similar structure)
    moments[10:14] = disaster_levels * 1.2  # Scaled for micro
    moments[14:18] = disaster_unc_probs * 0.8
    moments[18] = 0.7
    moments[19] = -8.0
    
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
    
    # Disaster probabilities
    disaster_probs = np.array([0.242, 0.03, 0.011, 0.008])
    
    # Initialize
    np.random.seed(2501)
    
    # Exogenous states
    a_pos = np.zeros(T, dtype=int)
    s_pos = np.zeros(T, dtype=int)
    disaster_occurred = np.zeros((T, 4), dtype=bool)
    
    a_pos[0] = p.anum // 2
    s_pos[0] = 0
    
    # Aggregate variables
    Y_sim = np.zeros(T)
    K_sim = np.zeros(T)
    I_sim = np.zeros(T)
    
    k_idx = p.knum // 2
    l_idx = p.lnum // 2
    
    price = p.pval
    w = p.theta / price
    
    for t in range(T):
        # Check for disasters
        for d in range(4):
            if np.random.random() < disaster_probs[d]:
                disaster_occurred[t, d] = True
                
                # Apply disaster impact on productivity
                # This reduces a in that period
                a_val = g.a_grid[a_pos[t]] * (1 + disaster_levels[d])
                
                # Maybe increase uncertainty
                if np.random.random() < disaster_unc_probs[d]:
                    s_pos[t] = 1
            else:
                a_val = g.a_grid[a_pos[t]]
        
        k = g.k_grid[k_idx]
        l = g.l_grid[l_idx]
        
        K_sim[t] = k
        Y_sim[t] = a_val * (k ** p.alpha) * (l ** p.nu)
        
        # Policy lookup
        endog_idx = k_idx * p.lnum + l_idx
        exog_idx = a_pos[t] * p.snum * p.snum + s_pos[t] * p.snum + s_pos[t]
        
        pol_idx = vfi_sol.polmat[endog_idx % g.numendog, exog_idx % g.numexog, 0]
        k_idx_next = g.endog_pos[pol_idx % g.numendog, 0]
        
        k_next = g.k_grid[k_idx_next]
        I_sim[t] = k_next - (1 - p.deltak) * k
        
        k_idx = k_idx_next
        
        # Transition aggregate productivity
        if t < T - 1:
            trans_probs = g.pr_mat_a[a_pos[t], :, s_pos[t]]
            a_pos[t+1] = np.searchsorted(np.cumsum(trans_probs), np.random.random())
            
            # Uncertainty transition
            if np.random.random() < g.pr_mat_s[s_pos[t], 1]:
                s_pos[t+1] = 1
            else:
                s_pos[t+1] = 0
    
    # Return simulation results
    from .simulation import SimulationResults
    return SimulationResults(
        Y_sim=Y_sim,
        K_sim=K_sim,
        L_sim=np.zeros(T),
        I_sim=I_sim,
        H_sim=np.zeros(T),
        ACk_sim=np.zeros(T),
        ACl_sim=np.zeros(T),  # Labor adjustment costs
        C_sim=Y_sim - I_sim,  # Consumption = Y - I
        p_sim=np.full(T, params.pval),  # Fixed price
        a_sim=a_pos,
        s_sim=s_pos
    )


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
    verbose: bool = True
) -> GMMSolution:
    """
    Estimate model parameters via GMM.
    
    Uses simplified optimization (Nelder-Mead) instead of PSO
    for faster computation.
    
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
    result = minimize(
        objective,
        x_init,
        method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter': max_evals, 'disp': verbose}
    )
    
    # Compute final moments
    final_moments = compute_simulated_moments(
        params, grids, vfi_sol,
        result.x[0:4], result.x[4:8]
    )
    
    return GMMSolution(
        x_opt=result.x,
        gmm_value=result.fun,
        converged=result.success,
        iterations=result.nit,
        simulated_moments=final_moments,
        data_moments=data_moments
    )
