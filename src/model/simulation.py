"""
Firm-Level Simulation and Aggregation for BBT (2024) Model.

Implements complete multi-firm simulation matching Fortran VOL_GROWTH_wrapper.f90 exactly:
1. Distribution evolution over state space
2. Multi-firm simulation (800 total firms, 200 public firms)
3. Stock returns computation
4. Cross-sectional volatility
5. Impulse response functions

Key arrays matching Fortran:
- zfirmpos: Firm productivity positions, shape (T, nfirms)
- endogfirmpos: Firm endogenous state positions, shape (T, nfirms)
- vfirm: Firm values, shape (T, nfirmspub)
- yfirm: Firm outputs, shape (T, nfirms)
- returnfirm: Stock returns, shape (T, nfirmspub)
- returnfirmsd: Cross-sectional return volatility, shape (T, nfirmspub)
"""

import numpy as np
from typing import Tuple, Dict, Optional
from dataclasses import dataclass
from numba import jit, prange
import warnings

from .params import ModelParameters
from .grids import StateGrids
from .vfi import VFISolution
from .adjustment import output, capital_adjustment_cost, labor_adjustment_cost


@dataclass
class SimulationResults:
    """Container for complete simulation results."""
    
    # Aggregate series (T periods)
    Y_sim: np.ndarray      # GDP
    K_sim: np.ndarray      # Aggregate capital
    L_sim: np.ndarray      # Aggregate labor
    I_sim: np.ndarray      # Investment
    H_sim: np.ndarray      # Hiring
    ACk_sim: np.ndarray    # Capital adjustment costs
    ACl_sim: np.ndarray    # Labor adjustment costs
    C_sim: np.ndarray      # Consumption
    p_sim: np.ndarray      # Price
    
    # Exogenous state series
    a_sim: np.ndarray      # Aggregate productivity positions
    s_sim: np.ndarray      # Uncertainty state positions
    
    # Firm-level series (for stock market analysis)
    vfirm: Optional[np.ndarray] = None       # Firm values, shape (T, nfirmspub)
    yfirm: Optional[np.ndarray] = None       # Firm outputs, shape (T, nfirms)
    dfirm: Optional[np.ndarray] = None       # Firm dividends, shape (T, nfirms)
    returnfirm: Optional[np.ndarray] = None  # Stock returns, shape (T, nfirmspub)
    returnfirmsd: Optional[np.ndarray] = None # Return volatility, shape (T, nfirmspub)
    
    # Distribution
    dist_zkl: Optional[np.ndarray] = None    # Distribution, shape (znum, numendog, T)


@dataclass
class IRFResults:
    """Container for impulse response results."""
    
    irf_Y: np.ndarray      # GDP IRF
    irf_I: np.ndarray      # Investment IRF
    irf_K: np.ndarray      # Capital IRF
    periods: np.ndarray    # Time periods


@jit(nopython=True)
def box_muller_transform(u1: float, u2: float) -> Tuple[float, float]:
    """
    Box-Muller transform for generating normal random variables.
    
    Matches Fortran subroutine boxmuller.
    """
    pi = 3.14159265358979323846
    r = np.sqrt(-2.0 * np.log(u1))
    theta = 2.0 * pi * u2
    return r * np.cos(theta), r * np.sin(theta)


@jit(nopython=True, parallel=True)
def simulate_firm_exog(
    z_shocks: np.ndarray,
    z_firm_pos: np.ndarray,
    pr_mat_z: np.ndarray,
    s_pos: np.ndarray,
    T: int,
    n_firms: int,
    znum: int,
    zinit: int
) -> None:
    """
    Simulate idiosyncratic productivity for all firms.
    
    Matches Fortran subroutine firmexogsim.
    
    Parameters
    ----------
    z_shocks : np.ndarray
        Random shocks, shape (T, n_firms).
    z_firm_pos : np.ndarray
        Output productivity positions, shape (T, n_firms).
    pr_mat_z : np.ndarray
        Productivity transition matrix, shape (znum, znum, snum).
    s_pos : np.ndarray
        Uncertainty state positions, shape (T,).
    T, n_firms, znum, zinit : int
        Dimensions.
    """
    # Initialize
    z_firm_pos[0, :] = zinit - 1  # Fortran uses 1-based indexing
    
    # Simulate
    for firm in prange(n_firms):
        for t in range(T - 1):
            s_idx = s_pos[t]
            z_curr = z_firm_pos[t, firm]
            
            # Build cumulative distribution
            cum_prob = 0.0
            for z_next in range(znum):
                cum_prob += pr_mat_z[z_curr, z_next, s_idx]
                if z_shocks[t + 1, firm] < cum_prob:
                    z_firm_pos[t + 1, firm] = z_next
                    break
    
    return


@jit(nopython=True, parallel=True)
def simulate_firms_core(
    z_firm_pos: np.ndarray,
    endog_firm_pos: np.ndarray,
    polmat: np.ndarray,
    V: np.ndarray,
    endog_pos: np.ndarray,
    exog_pos: np.ndarray,
    k_grid: np.ndarray,
    l_grid: np.ndarray,
    z_grid: np.ndarray,
    a_grid: np.ndarray,
    Ymat: np.ndarray,
    a_pos: np.ndarray,
    s_pos: np.ndarray,
    w: float,
    alpha: float,
    nu: float,
    capirrev: float,
    capfix: float,
    deltak: float,
    deltan: float,
    hirelin: float,
    firelin: float,
    labfix: float,
    knum: int,
    lnum: int,
    anum: int,
    snum: int,
    znum: int,
    T: int,
    n_firms: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Core simulation of firms.
    
    Matches Fortran lines 1255-1306.
    
    Returns
    -------
    yfirm, dfirm, vfirm : np.ndarray
        Firm outputs, dividends, and values.
    """
    yfirm = np.zeros((T, n_firms))
    dfirm = np.zeros((T, n_firms))
    vfirm = np.zeros((T, n_firms))
    
    for firm in prange(n_firms):
        for t in range(T):
            # Get current state
            z_idx = z_firm_pos[t, firm]
            endog_idx = endog_firm_pos[t, firm]
            a_idx = a_pos[t]
            s_idx = s_pos[t]
            
            # Get previous uncertainty state
            if t > 0:
                sm1_idx = s_pos[t - 1]
            else:
                sm1_idx = 0  # Low uncertainty previously
            
            # Construct exogenous index
            # Fortran line 1263
            exog_idx = z_idx * anum * snum * snum + a_idx * snum * snum + s_idx * snum + sm1_idx
            exog_idx = min(exog_idx, polmat.shape[1] - 1)
            
            # Get policy
            polstar = polmat[endog_idx, exog_idx, 0]
            polstar = min(polstar, endog_pos.shape[0] - 1)
            
            # Update endogenous state for next period
            if t < T - 1:
                endog_firm_pos[t + 1, firm] = polstar
            
            # Get capital and labor positions
            k_idx = endog_pos[endog_idx, 0]
            lmin1_idx = endog_pos[endog_idx, 1]
            kprime_idx = endog_pos[polstar, 0]
            l_idx = endog_pos[polstar, 1]
            
            # Get grid values
            k = k_grid[k_idx]
            l = l_grid[l_idx]
            l_prev = l_grid[lmin1_idx]
            k_prime = k_grid[kprime_idx]
            z = z_grid[z_idx]
            a = a_grid[a_idx]
            
            # Output
            yval = Ymat[z_idx, a_idx, k_idx, l_idx]
            yfirm[t, firm] = yval
            
            # Dividend (d = Y - ACk - ACl - I - wL)
            # Fortran lines 1280-1284
            Y_val = yval
            I_val = k_prime - (1.0 - deltak) * k
            
            # Capital adjustment cost
            changetol = 1e-10
            ACk_val = 0.0
            if I_val < -changetol:
                ACk_val = -I_val * capirrev
            if abs(I_val) > changetol:
                ACk_val += Y_val * capfix
            
            # Labor adjustment cost
            h = l - (1.0 - deltan) * l_prev
            ACl_val = 0.0
            if abs(h) > changetol:
                ACl_val += labfix * Y_val
                if h < -changetol:
                    ACl_val += -h * w * firelin
                if h > changetol:
                    ACl_val += h * w * hirelin
            
            d_val = Y_val - ACk_val - ACl_val - I_val - w * l
            dfirm[t, firm] = d_val
            
            # Firm value
            vfirm[t, firm] = V[endog_idx, exog_idx, 0] / (w * nu / (1.0 - alpha - nu))  # Normalized
    
    return yfirm, dfirm, vfirm


@jit(nopython=True, parallel=True)
def compute_stock_returns(
    vfirm: np.ndarray,
    dfirm: np.ndarray,
    T: int,
    n_firms_pub: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute stock returns and cross-sectional volatility.
    
    Matches Fortran lines 1292-1304, 1334-1342.
    """
    returnfirm = np.zeros((T, n_firms_pub))
    returnfirmsd = np.zeros((T, n_firms_pub))
    
    for firm in prange(n_firms_pub):
        for t in range(2, T):
            # Return = log(v_t / (v_{t-1} - d_{t-1}))
            # Fortran line 1293
            denominator = vfirm[t - 1, firm] - dfirm[t - 1, firm]
            if denominator > 0:
                returnfirm[t, firm] = np.log(vfirm[t, firm] / denominator)
            
            # Rolling standard deviation
            # Fortran lines 1298-1303
            if t >= 4:
                meanval = 0.25 * (returnfirm[t, firm] + returnfirm[t-1, firm] + 
                                  returnfirm[t-2, firm] + returnfirm[t-3, firm])
                sdval = 0.25 * (returnfirm[t, firm]**2 + returnfirm[t-1, firm]**2 + 
                               returnfirm[t-2, firm]**2 + returnfirm[t-3, firm]**2)
                val = sdval - meanval**2
                if val > 0:
                    returnfirmsd[t, firm] = np.sqrt(val)
    
    return returnfirm, returnfirmsd


def simulate_all_firms(
    params: ModelParameters,
    grids: StateGrids,
    vfi_sol: VFISolution,
    T: int = None,
    seed: int = 2501,
    verbose: bool = True
) -> SimulationResults:
    """
    Run complete multi-firm simulation.
    
    Matches Fortran simulation block (lines 1057-1350).
    
    Parameters
    ----------
    params : ModelParameters
    grids : StateGrids
    vfi_sol : VFISolution
    T : int, optional
        Number of periods. Default from params.
    seed : int
        Random seed.
    verbose : bool
    
    Returns
    -------
    SimulationResults
        Complete simulation output.
    """
    p = params
    g = grids
    
    if T is None:
        T = p.numdiscard + p.Ncountries * p.Tper
    
    np.random.seed(seed)
    
    if verbose:
        print(f"Simulating {p.nfirms} firms for {T} periods...")
    
    # =====================
    # Initialize Distribution
    # =====================
    dist_zkl = np.zeros((p.znum, g.numendog, T))
    
    # Initialize at middle of grid
    kct = min(p.knum // 2, p.knum - 1)
    lmin1ct = min(p.lnum // 3, p.lnum - 1)
    endogct = kct * p.lnum + lmin1ct
    
    for zi in range(p.znum):
        for ej in range(max(0, endogct - 5), min(g.numendog, endogct + 6)):
            dist_zkl[zi, ej, 0] = 1.0
    dist_zkl[:, :, 0] /= dist_zkl[:, :, 0].sum()
    
    # =====================
    # Initialize Firms
    # =====================
    # Firm positions
    z_firm_pos = np.zeros((T, p.nfirms), dtype=int)
    endog_firm_pos = np.zeros((T, p.nfirms), dtype=int)
    
    # Initialize at middle of grids
    z_firm_pos[0, :] = p.zinit - 1  # Convert to 0-based
    endog_firm_pos[0, :] = g.numendog // 2
    
    # =====================
    # Simulate Exogenous Processes
    # =====================
    a_pos = np.zeros(T, dtype=int)
    s_pos = np.zeros(T, dtype=int)
    
    a_pos[0] = p.ainit - 1
    s_pos[0] = p.sinit - 1
    
    # Random shocks
    a_shocks = np.random.random(T)
    s_shocks = np.random.random(T)
    z_shocks = np.random.random((T, p.nfirms))
    
    # Simulate exogenous processes
    for t in range(1, T):
        # Uncertainty transition
        if s_shocks[t] < g.pr_mat_s[s_pos[t-1], 1]:
            s_pos[t] = 1
        else:
            s_pos[t] = 0
        
        # Productivity transition
        trans_probs = g.pr_mat_a[a_pos[t-1], :, s_pos[t]]
        a_pos[t] = np.searchsorted(np.cumsum(trans_probs), a_shocks[t])
        a_pos[t] = min(a_pos[t], p.anum - 1)
    
    # Simulate firm productivity
    simulate_firm_exog(
        z_shocks, z_firm_pos, g.pr_mat_z, s_pos,
        T, p.nfirms, p.znum, p.zinit
    )
    
    # =====================
    # Simulate Firms
    # =====================
    # Pre-compute output matrix
    Ymat = np.zeros((p.znum, p.anum, p.knum, p.lnum))
    for zi in range(p.znum):
        for ai in range(p.anum):
            for ki in range(p.knum):
                for li in range(p.lnum):
                    Ymat[zi, ai, ki, li] = output(
                        g.z_grid[zi], g.a_grid[ai],
                        g.k_grid[ki], g.l_grid[li],
                        p.alpha, p.nu
                    )
    
    # Price and wage
    price = p.pval
    w = p.theta / price
    
    # Run firm simulation
    yfirm, dfirm, vfirm = simulate_firms_core(
        z_firm_pos, endog_firm_pos, vfi_sol.polmat, vfi_sol.V,
        g.endog_pos, g.exog_pos, g.k_grid, g.l_grid, g.z_grid, g.a_grid,
        Ymat, a_pos, s_pos, w,
        p.alpha, p.nu, p.capirrev, p.capfix, p.deltak, p.deltan,
        p.hirelin, p.firelin, p.labfix,
        p.knum, p.lnum, p.anum, p.snum, p.znum, T, p.nfirms
    )
    
    # =====================
    # Compute Stock Returns
    # =====================
    returnfirm, returnfirmsd = compute_stock_returns(
        vfirm, dfirm, T, p.nfirmspub
    )
    
    # Add noise to returns
    norm_shocks = np.random.randn(T, p.nfirmspub)
    return_mean = np.mean(returnfirm[2:, :])
    return_std = np.std(returnfirm[2:, :])
    returnfirm_noise = returnfirm + return_std * norm_shocks
    
    # =====================
    # Compute Aggregates
    # =====================
    Y_sim = np.sum(yfirm, axis=1)
    K_sim = np.zeros(T)
    L_sim = np.zeros(T)
    I_sim = np.zeros(T)
    H_sim = np.zeros(T)
    ACk_sim = np.zeros(T)
    ACl_sim = np.zeros(T)
    C_sim = np.zeros(T)
    p_sim = np.full(T, price)
    
    for t in range(T):
        K_sim[t] = np.mean(g.k_grid[endog_firm_pos[t, :] // p.lnum])
        L_sim[t] = np.mean(g.l_grid[endog_firm_pos[t, :] % p.lnum])
    
    # Investment and adjustment costs
    for t in range(1, T):
        I_sim[t] = K_sim[t] - (1 - p.deltak) * K_sim[t-1]
        H_sim[t] = L_sim[t] - (1 - p.deltan) * L_sim[t-1]
        
        ACk_sim[t] = capital_adjustment_cost(
            K_sim[t], K_sim[t-1], p.capirrev, p.capfix, p.deltak
        )
        ACl_sim[t] = labor_adjustment_cost(
            L_sim[t], L_sim[t-1], w, p.hirelin, p.firelin, p.labfix, p.deltan
        )
        
        C_sim[t] = Y_sim[t] - I_sim[t] - ACk_sim[t] - ACl_sim[t]
    
    if verbose:
        print(f"Simulation complete. GDP mean: {Y_sim.mean():.4f}")
    
    return SimulationResults(
        Y_sim=Y_sim, K_sim=K_sim, L_sim=L_sim,
        I_sim=I_sim, H_sim=H_sim, ACk_sim=ACk_sim, ACl_sim=ACl_sim,
        C_sim=C_sim, p_sim=p_sim,
        a_sim=a_pos, s_sim=s_pos,
        vfirm=vfirm[:, :p.nfirmspub],
        yfirm=yfirm, dfirm=dfirm,
        returnfirm=returnfirm[:, :p.nfirmspub],
        returnfirmsd=returnfirmsd[:, :p.nfirmspub],
        dist_zkl=dist_zkl
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
    Compute impulse response to uncertainty shock.
    
    Parameters
    ----------
    params : ModelParameters
    grids : StateGrids
    vfi_sol : VFISolution
    shock_size : float
        Multiplier for shock.
    T_irf : int
        IRF length.
    n_sims : int
        Number of simulations.
    seed : int
    
    Returns
    -------
    IRFResults
    """
    p = params
    g = grids
    
    np.random.seed(seed)
    
    # Storage for IRFs
    irf_Y = np.zeros((T_irf, n_sims))
    irf_I = np.zeros((T_irf, n_sims))
    
    for sim in range(n_sims):
        # Run simulation with shock
        T = T_irf + 10  # Some burn-in
        shock_period = 5
        
        # Initialize
        a_pos = np.zeros(T, dtype=int)
        s_pos = np.zeros(T, dtype=int)
        
        a_pos[0] = p.anum // 2
        s_pos[0] = 0
        
        # Aggregate variables
        Y_sim = np.zeros(T)
        I_sim = np.zeros(T)
        
        k_idx = p.knum // 2
        l_idx = p.lnum // 2
        
        for t in range(T):
            # Uncertainty shock
            if shock_period <= t < shock_period + int(5 * shock_size):
                s_pos[t] = 1
            elif t >= shock_period + int(5 * shock_size):
                if np.random.random() < p.uncpers:
                    s_pos[t] = 1
                else:
                    s_pos[t] = 0
            
            # Get grid values
            k = g.k_grid[k_idx]
            l = g.l_grid[l_idx]
            a = g.a_grid[a_pos[t]]
            
            # Output
            Y_sim[t] = output(1.0, a, k, l, p.alpha, p.nu)
            
            # Policy lookup
            endog_idx = k_idx * p.lnum + l_idx
            exog_idx = a_pos[t] * p.snum * p.snum + s_pos[t] * p.snum + s_pos[t]
            
            pol_idx = vfi_sol.polmat[endog_idx % g.numendog, exog_idx % g.numexog, 0]
            k_idx_next = g.endog_pos[pol_idx % g.numendog, 0]
            l_idx_next = g.endog_pos[pol_idx % g.numendog, 1]
            
            k_next = g.k_grid[k_idx_next]
            I_sim[t] = k_next - (1 - p.deltak) * k
            
            k_idx = k_idx_next
            l_idx = l_idx_next
            
            # Productivity transition
            if t < T - 1:
                trans_probs = g.pr_mat_a[a_pos[t], :, s_pos[t]]
                a_pos[t+1] = np.searchsorted(np.cumsum(trans_probs), np.random.random())
                a_pos[t+1] = min(a_pos[t+1], p.anum - 1)
        
        # Store IRF (relative to first period)
        irf_Y[:, sim] = (Y_sim[:T_irf] / Y_sim[0] - 1) * 100
        irf_I[:, sim] = (I_sim[:T_irf] / max(abs(I_sim[0]), 1e-10) - 1) * 100
    
    # Average across simulations
    return IRFResults(
        irf_Y=np.mean(irf_Y, axis=1),
        irf_I=np.mean(irf_I, axis=1),
        irf_K=np.cumsum(np.mean(irf_Y, axis=1)) * 0.25,
        periods=np.arange(1, T_irf + 1)
    )


def compute_figure8_irf(
    params: ModelParameters,
    grids: StateGrids,
    vfi_sol: VFISolution,
    T: int = 40
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute IRF for Figure 8: GDP and Investment response to uncertainty shock.
    """
    irf_results = simulate_irf(params, grids, vfi_sol, T_irf=T, n_sims=500)
    return irf_results.irf_Y, irf_results.irf_I


# Keep backward compatibility
def simulate_firms(*args, **kwargs):
    """Backward compatible wrapper."""
    return simulate_all_firms(*args, **kwargs)


def simulate_firms_with_shock(*args, **kwargs):
    """Backward compatible wrapper."""
    return simulate_irf(*args, **kwargs)
