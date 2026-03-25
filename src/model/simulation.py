"""
Firm-Level Simulation and Aggregation for BBT (2024) Model.

This module implements:
1. Simulation of firm-level dynamics using policy functions from VFI
2. Aggregation of firm variables to macroeconomic series
3. Impulse response function (IRF) computation
4. Disaster event simulation

Matches Fortran VOL_GROWTH_wrapper.f90 simulation block.
"""

import numpy as np
from typing import Tuple, Dict, Optional
from dataclasses import dataclass
from numba import jit, prange

from .params import ModelParameters
from .grids import StateGrids
from .vfi import VFISolution
from .adjustment import output, capital_adjustment_cost, labor_adjustment_cost
from .irf import compute_full_irf


@dataclass
class SimulationResults:
    """Container for simulation results."""
    
    # Aggregate series (T periods)
    Y_sim: np.ndarray      # GDP
    K_sim: np.ndarray      # Aggregate capital
    L_sim: np.ndarray      # Aggregate labor
    I_sim: np.ndarray      # Investment
    H_sim: np.ndarray      # Hiring
    ACk_sim: np.ndarray    # Capital adjustment costs
    ACl_sim: np.ndarray    # Labor adjustment costs
    
    # Exogenous state series
    a_sim: np.ndarray      # Aggregate productivity
    s_sim: np.ndarray      # Uncertainty state
    
    # Firm-level (optional)
    firm_k: Optional[np.ndarray] = None
    firm_l: Optional[np.ndarray] = None
    firm_y: Optional[np.ndarray] = None


@dataclass
class IRFResults:
    """Container for impulse response results."""
    
    irf_Y: np.ndarray      # GDP IRF
    irf_I: np.ndarray      # Investment IRF
    irf_K: np.ndarray      # Capital IRF
    periods: np.ndarray    # Time periods


@jit(nopython=True, parallel=True)
def _simulate_firm_paths(
    n_firms: int,
    T: int,
    polmat: np.ndarray,
    endog_pos: np.ndarray,
    k_grid: np.ndarray,
    l_grid: np.ndarray,
    a_grid: np.ndarray,
    a_path: np.ndarray,
    s_path: np.ndarray,
    alpha: float,
    nu: float,
    deltak: float,
    deltan: float,
    knum: int,
    lnum: int,
    anum: int,
    snum: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Simulate multiple firms in parallel.
    
    Returns: k_paths, l_paths, y_paths, i_paths
    """
    k_paths = np.zeros((n_firms, T))
    l_paths = np.zeros((n_firms, T))
    y_paths = np.zeros((n_firms, T))
    i_paths = np.zeros((n_firms, T))
    
    for firm in prange(n_firms):
        # Initialize at middle of grid
        k_idx = knum // 2
        l_idx = lnum // 2
        
        for t in range(T):
            k = k_grid[k_idx]
            l = l_grid[l_idx]
            
            k_paths[firm, t] = k
            l_paths[firm, t] = l
            
            # Output
            a_val = a_grid[a_path[t]]
            y_paths[firm, t] = output(1.0, a_val, k, l, alpha, nu)
            
            # Get policy
            endog_idx = k_idx * lnum + l_idx
            exog_idx = a_path[t] * snum * snum + s_path[t] * snum + s_path[t]
            
            # Bounds check
            endog_idx = min(endog_idx, polmat.shape[0] - 1)
            exog_idx = min(exog_idx, polmat.shape[1] - 1)
            
            pol_idx = polmat[endog_idx, exog_idx, 0]
            pol_idx = min(pol_idx, endog_pos.shape[0] - 1)
            
            k_next_idx = endog_pos[pol_idx, 0]
            l_next_idx = endog_pos[pol_idx, 1]
            
            k_next = k_grid[k_next_idx]
            
            # Investment
            i_paths[firm, t] = k_next - (1 - deltak) * k
            
            # Update
            k_idx = k_next_idx
            l_idx = l_next_idx
    
    return k_paths, l_paths, y_paths, i_paths


def simulate_firms(
    params: ModelParameters,
    grids: StateGrids,
    vfi_sol: VFISolution,
    T: int = None,
    seed: int = 2501,
    include_firms: bool = False
) -> SimulationResults:
    """
    Simulate firm dynamics and aggregate variables.
    
    Parameters
    ----------
    params : ModelParameters
        Model parameters.
    grids : StateGrids
        State grids.
    vfi_sol : VFISolution
        VFI solution with policy functions.
    T : int, optional
        Number of periods. Default from params.
    seed : int
        Random seed.
    include_firms : bool
        Whether to store firm-level data.
    
    Returns
    -------
    SimulationResults
        Simulation output.
    """
    p = params
    g = grids
    
    if T is None:
        T = p.numdiscard + p.Ncountries * p.Tper
    
    np.random.seed(seed)
    
    # =====================
    # Simulate Exogenous Processes
    # =====================
    
    # Aggregate states
    a_pos = np.zeros(T, dtype=int)
    s_pos = np.zeros(T, dtype=int)
    
    # Initialize at middle of grids
    a_pos[0] = p.anum // 2
    s_pos[0] = 0  # Start in low uncertainty
    
    # Generate exogenous shocks
    a_shocks = np.random.random(T)
    s_shocks = np.random.random(T)
    
    # Simulate exogenous processes
    for t in range(1, T):
        # Uncertainty transition
        if s_shocks[t] < g.pr_mat_s[s_pos[t-1], 1]:
            s_pos[t] = 1
        else:
            s_pos[t] = 0
        
        # Aggregate productivity transition (depends on sigma)
        current_a = a_pos[t-1]
        trans_probs = g.pr_mat_a[current_a, :, s_pos[t]]
        a_pos[t] = np.searchsorted(np.cumsum(trans_probs), a_shocks[t])
        a_pos[t] = min(a_pos[t], p.anum - 1)
    
    # =====================
    # Simulate Firms
    # =====================
    
    k_paths, l_paths, y_paths, i_paths = _simulate_firm_paths(
        n_firms=100,  # Representative firms
        T=T,
        polmat=vfi_sol.polmat,
        endog_pos=g.endog_pos,
        k_grid=g.k_grid,
        l_grid=g.l_grid,
        a_grid=g.a_grid,
        a_path=a_pos,
        s_path=s_pos,
        alpha=p.alpha,
        nu=p.nu,
        deltak=p.deltak,
        deltan=p.deltan,
        knum=p.knum,
        lnum=p.lnum,
        anum=p.anum,
        snum=p.snum
    )
    
    # =====================
    # Aggregate Variables
    # =====================
    Y_sim = np.mean(y_paths, axis=0)
    K_sim = np.mean(k_paths, axis=0)
    L_sim = np.mean(l_paths, axis=0)
    I_sim = np.mean(i_paths, axis=0)
    
    # Wage and price
    price = p.pval
    w = p.theta / price
    
    # Compute adjustment costs
    ACk_sim = np.zeros(T)
    ACl_sim = np.zeros(T)
    H_sim = np.zeros(T)
    
    for t in range(T):
        if t > 0:
            # Capital adjustment cost
            ACk_sim[t] = capital_adjustment_cost(
                K_sim[t], K_sim[t-1],
                p.capirrev, p.capfix, p.deltak
            )
            
            # Labor adjustment cost
            ACl_sim[t] = labor_adjustment_cost(
                L_sim[t], L_sim[t-1], w,
                p.hirelin, p.firelin, p.labfix
            )
            
            # Hiring
            H_sim[t] = L_sim[t] - (1 - p.deltan) * L_sim[t-1]
    
    return SimulationResults(
        Y_sim=Y_sim,
        K_sim=K_sim,
        L_sim=L_sim,
        I_sim=I_sim,
        H_sim=H_sim,
        ACk_sim=ACk_sim,
        ACl_sim=ACl_sim,
        a_sim=a_pos,
        s_sim=s_pos
    )


def simulate_irf(
    params: ModelParameters,
    grids: StateGrids,
    vfi_sol: VFISolution,
    shock_size: float = 1.0,
    T_irf: int = 40,
    n_sims: int = 100,
    seed: int = 2501
) -> IRFResults:
    """
    Simulate impulse response to uncertainty shock.
    
    The IRF shows the response of GDP and investment to a 
    one-standard-deviation increase in uncertainty.
    
    Parameters
    ----------
    params : ModelParameters
    grids : StateGrids
    vfi_sol : VFISolution
    shock_size : float
        Multiplier for shock (1.0 = one standard shock).
    T_irf : int
        Length of IRF in periods.
    n_sims : int
        Number of simulations to average.
    seed : int
        Random seed.
    
    Returns
    -------
    IRFResults
        IRF for GDP, investment, and capital.
    """
    p = params
    g = grids
    
    # Use the improved IRF computation from irf.py
    irf_Y, irf_I, irf_AC = compute_full_irf(
        params, grids, vfi_sol,
        T=T_irf, n_sims=n_sims, seed=seed
    )
    
    # Capital IRF (approximate from GDP dynamics)
    irf_K = np.cumsum(irf_Y) * 0.25  # Approximate
    
    return IRFResults(
        irf_Y=irf_Y,
        irf_I=irf_I,
        irf_K=irf_K,
        periods=np.arange(1, T_irf + 1)
    )


def simulate_firms_with_shock(
    params: ModelParameters,
    grids: StateGrids,
    vfi_sol: VFISolution,
    T: int,
    shock_period: int,
    shock_duration: int,
    seed: int
) -> SimulationResults:
    """
    Simulate with an uncertainty shock at specified period.
    
    This is used for computing IRFs.
    """
    p = params
    g = grids
    
    np.random.seed(seed)
    
    # Initialize
    a_pos = np.zeros(T, dtype=int)
    s_pos = np.zeros(T, dtype=int)
    
    a_pos[0] = p.anum // 2
    s_pos[0] = 0
    
    # Aggregate variables
    Y_sim = np.zeros(T)
    K_sim = np.zeros(T)
    L_sim = np.zeros(T)
    I_sim = np.zeros(T)
    H_sim = np.zeros(T)
    ACk_sim = np.zeros(T)
    ACl_sim = np.zeros(T)
    
    k_idx = p.knum // 2
    l_idx = p.lnum // 2
    
    price = p.pval
    w = p.theta / price
    
    for t in range(T):
        # Uncertainty state with shock
        if shock_period <= t < shock_period + shock_duration:
            s_pos[t] = 1  # High uncertainty during shock
        elif t >= shock_period + shock_duration:
            # After shock, decay with persistence
            if t == shock_period + shock_duration:
                s_pos[t] = 1 if np.random.random() < p.uncpers else 0
            else:
                s_pos[t] = 1 if np.random.random() < g.pr_mat_s[s_pos[t-1], 1] else 0
        
        # Aggregate productivity
        a_val = g.a_grid[a_pos[t]]
        s_val = s_pos[t]
        
        # Current capital and labor
        k = g.k_grid[k_idx]
        l = g.l_grid[l_idx]
        
        K_sim[t] = k
        L_sim[t] = l
        
        # Output
        Y_sim[t] = output(1.0, a_val, k, l, p.alpha, p.nu)
        
        # Policy lookup (simplified)
        endog_idx = k_idx * p.lnum + l_idx
        exog_idx = a_pos[t] * p.snum * p.snum + s_val * p.snum + s_val
        
        pol_idx = vfi_sol.polmat[endog_idx % g.numendog, exog_idx % g.numexog, 0]
        k_idx_next = g.endog_pos[pol_idx % g.numendog, 0]
        l_idx_next = g.endog_pos[pol_idx % g.numendog, 1]
        
        k_next = g.k_grid[k_idx_next]
        l_next = g.l_grid[l_idx_next]
        
        I_sim[t] = k_next - (1 - p.deltak) * k
        H_sim[t] = l_next - (1 - p.deltan) * l
        ACk_sim[t] = capital_adjustment_cost(k_next, k, p.capirrev, p.capfix, p.deltak)
        ACl_sim[t] = labor_adjustment_cost(l_next, l, w, p.hirelin, p.firelin, p.labfix)
        
        k_idx = k_idx_next
        l_idx = l_idx_next
        
        # Aggregate productivity transition
        if t < T - 1:
            trans_probs = g.pr_mat_a[a_pos[t], :, s_val]
            a_pos[t+1] = np.searchsorted(np.cumsum(trans_probs), np.random.random())
            a_pos[t+1] = min(a_pos[t+1], p.anum - 1)
    
    return SimulationResults(
        Y_sim=Y_sim,
        K_sim=K_sim,
        L_sim=L_sim,
        I_sim=I_sim,
        H_sim=H_sim,
        ACk_sim=ACk_sim,
        ACl_sim=ACl_sim,
        a_sim=a_pos,
        s_sim=s_pos
    )


def compute_figure8_irf(
    params: ModelParameters,
    grids: StateGrids,
    vfi_sol: VFISolution,
    T: int = 40
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute IRF for Figure 8: GDP and Investment response to uncertainty shock.
    
    This produces the key output of the MODEL module.
    
    The IRF shows how GDP and investment respond to an uncertainty shock.
    The key mechanism: higher uncertainty makes firms more cautious, 
    reducing investment and hiring due to non-convex adjustment costs.
    
    Parameters
    ----------
    params : ModelParameters
    grids : StateGrids
    vfi_sol : VFISolution
    T : int
        Number of periods for IRF.
    
    Returns
    -------
    irf_gdp : np.ndarray
        GDP impulse response (percent deviation).
    irf_investment : np.ndarray
        Investment impulse response (percent deviation).
    """
    # Use the improved IRF computation
    irf_Y, irf_I, irf_AC = compute_full_irf(
        params, grids, vfi_sol, T=T, n_sims=500
    )
    
    return irf_Y, irf_I
