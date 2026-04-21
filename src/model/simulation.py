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

CRITICAL FIX: Use MT19937 random generator matching Fortran's random_number
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


def create_mt19937_rng(seed: int = 2501):
    """
    Create a random number generator matching Fortran's random_number.
    
    Fortran uses a Mersenne Twister by default in modern compilers.
    Python's RandomState uses MT19937.
    """
    rng = np.random.RandomState(seed)
    return rng


def draw_uniform_matrix_fortran(rng: np.random.RandomState, nrow: int, ncol: int) -> np.ndarray:
    """
    Draw U(0,1) matrix emulating Fortran column-major fill order.
    """
    flat = rng.random(nrow * ncol)
    return flat.reshape((nrow, ncol), order="F")


def box_muller_matrix_from_uniforms(u1: np.ndarray, u2: np.ndarray) -> np.ndarray:
    """
    Generate normal shocks via Box-Muller from two U(0,1) matrices.
    """
    u1c = np.clip(u1, 1e-12, 1.0 - 1e-12)
    r = np.sqrt(-2.0 * np.log(u1c))
    theta = 2.0 * np.pi * u2
    return r * np.cos(theta)


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
    returnfirm_noise: Optional[np.ndarray] = None   # Noisy returns, shape (T, nfirmspub)
    returnfirmsd_noise: Optional[np.ndarray] = None # Noisy return vol, shape (T, nfirmspub)
    
    # Distribution
    dist_zkl: Optional[np.ndarray] = None    # Distribution, shape (znum, numendog, T)
    
    # Optional event indicators used by GMM/IV moment construction
    disaster_indicators: Optional[np.ndarray] = None  # shape (T, 4)


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
        for t in range(1, T):
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
    verbose: bool = True,
    a_pos_override: Optional[np.ndarray] = None,
    s_pos_override: Optional[np.ndarray] = None,
    z_shocks_override: Optional[np.ndarray] = None,
    norm_shocks_override: Optional[np.ndarray] = None,
    disaster_indicators: Optional[np.ndarray] = None,
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
        Random seed (Fortran uses 2501 for simulation).
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
    
    # CRITICAL FIX: Use MT19937 matching Fortran's random_number
    rng = create_mt19937_rng(seed)
    
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
    if a_pos_override is not None and s_pos_override is not None:
        a_pos = np.asarray(a_pos_override, dtype=int).copy()
        s_pos = np.asarray(s_pos_override, dtype=int).copy()
        if len(a_pos) != T or len(s_pos) != T:
            raise ValueError("Override exogenous paths must have length T")
        a_pos = np.clip(a_pos, 0, p.anum - 1)
        s_pos = np.clip(s_pos, 0, p.snum - 1)
    else:
        a_pos = np.zeros(T, dtype=int)
        s_pos = np.zeros(T, dtype=int)
        
        a_pos[0] = p.ainit - 1
        s_pos[0] = p.sinit - 1
        
        # Random shocks - CRITICAL FIX: use rng instead of np.random
        a_shocks = rng.random(T)
        s_shocks = rng.random(T)
        
        # Simulate exogenous processes
        for t in range(1, T):
            # Uncertainty transition
            cdf_s = np.cumsum(g.pr_mat_s[s_pos[t - 1], :])
            s_pos[t] = int(np.searchsorted(cdf_s, s_shocks[t]))
            s_pos[t] = min(s_pos[t], p.snum - 1)
            
            # Productivity transition
            # Fortran uses current-period uncertainty state from t-1 transition.
            trans_probs = g.pr_mat_a[a_pos[t - 1], :, s_pos[t - 1]]
            a_pos[t] = np.searchsorted(np.cumsum(trans_probs), a_shocks[t])
            a_pos[t] = min(a_pos[t], p.anum - 1)

    if z_shocks_override is not None:
        z_shocks = np.asarray(z_shocks_override, dtype=np.float64)
        if z_shocks.shape != (T, p.nfirms):
            raise ValueError("z_shocks_override must have shape (T, nfirms)")
    else:
        z_shocks = draw_uniform_matrix_fortran(rng, T, p.nfirms)
    
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
    
    # Add noise to returns - CRITICAL FIX: use rng
    if norm_shocks_override is not None:
        norm_shocks = np.asarray(norm_shocks_override, dtype=np.float64)
        if norm_shocks.shape != (T, p.nfirmspub):
            raise ValueError("norm_shocks_override must have shape (T, nfirmspub)")
    else:
        u1 = draw_uniform_matrix_fortran(rng, T, p.nfirmspub)
        u2 = draw_uniform_matrix_fortran(rng, T, p.nfirmspub)
        norm_shocks = box_muller_matrix_from_uniforms(u1, u2)
    # Fortran uses population moments over t=2..numper-1 (1-based),
    # i.e. idx 1..T-2 (0-based), with divisor ct (not ct-1).
    return_slice = returnfirm[1:T-1, :]
    return_mean = np.mean(return_slice)
    return_std = np.sqrt(np.mean(return_slice ** 2) - return_mean ** 2)
    if not np.isfinite(return_std) or return_std < 0:
        return_std = 0.0
    returnfirm_noise = returnfirm + return_std * norm_shocks
    returnfirmsd_noise = np.zeros_like(returnfirmsd)

    # Recompute within-firm rolling return volatility using noisy returns
    # (Fortran lines 1334-1340).
    # Fortran loop t=4..numper-1 (1-based) => idx 3..T-2 (0-based).
    for t in range(3, T - 1):
        wnd = returnfirm_noise[t - 3:t + 1, :]
        meanval = np.mean(wnd, axis=0)
        sdval = np.mean(wnd ** 2, axis=0)
        var = sdval - meanval ** 2
        var[var < 0] = 0.0
        returnfirmsd_noise[t, :] = np.sqrt(var)
    
    # =====================
    # Compute Aggregates (Fortran distribution-based accounting)
    # =====================
    Y_sim = np.zeros(T)
    K_sim = np.zeros(T)
    L_sim = np.zeros(T)
    I_sim = np.zeros(T)
    H_sim = np.zeros(T)
    ACk_sim = np.zeros(T)
    ACl_sim = np.zeros(T)
    C_sim = np.zeros(T)
    p_sim = np.full(T, price)
    K_sim[0] = (p.kbarmin + p.kbarmax) / 2.0

    # Fortran main simulation loop is t=1..numper-1 (1-based).
    # Here: t=0..T-2 (0-based); last period remains default zeros.
    for t in range(T - 1):
        act = int(a_pos[t])
        sct = int(s_pos[t])
        smin1ct = int(s_pos[t - 1]) if t > 0 else 0

        dist_t = dist_zkl[:, :, t]
        if t + 1 < T:
            dist_zkl[:, :, t + 1] = 0.0

        Y_val = 0.0
        I_val = 0.0
        ACk_val = 0.0
        ACl_val = 0.0
        Kprime_val = 0.0
        H_val = 0.0
        L_val = 0.0

        for zct in range(p.znum):
            for endogct in range(g.numendog):
                weight = dist_t[zct, endogct]
                if weight <= p.disttol:
                    continue

                exog_idx = (
                    zct * p.anum * p.snum * p.snum
                    + act * p.snum * p.snum
                    + sct * p.snum
                    + smin1ct
                )
                exog_idx = min(exog_idx, vfi_sol.polmat.shape[1] - 1)
                polstar = int(vfi_sol.polmat[endogct, exog_idx, 0])

                kct = int(g.endog_pos[endogct, 0])
                lmin1ct = int(g.endog_pos[endogct, 1])
                kprimect = int(g.endog_pos[polstar, 0])
                lct = int(g.endog_pos[polstar, 1])

                k = g.k_grid[kct]
                l = g.l_grid[lct]
                l_prev = g.l_grid[lmin1ct]
                k_prime = g.k_grid[kprimect]
                y_val = Ymat[zct, act, kct, lct]
                i_val = k_prime - (1.0 - p.deltak) * k
                ack_val = capital_adjustment_cost(k_prime, k, p.capirrev, p.capfix, p.deltak)
                acl_val = labor_adjustment_cost(
                    l, l_prev, w, p.hirelin, p.firelin, p.labfix, p.deltan
                )

                Y_val += weight * y_val
                I_val += weight * i_val
                ACk_val += weight * ack_val
                ACl_val += weight * acl_val
                Kprime_val += weight * k_prime
                H_val += weight * (l - (1.0 - p.deltan) * l_prev)
                L_val += weight * l

                if t + 1 < T:
                    for zprime in range(p.znum):
                        dist_zkl[zprime, polstar, t + 1] += g.pr_mat_z[zct, zprime, sct] * weight

        if t + 1 < T:
            denom = dist_zkl[:, :, t + 1].sum()
            if denom > 0:
                dist_zkl[:, :, t + 1] /= denom
            else:
                dist_zkl[:, :, t + 1] = dist_t

        Y_sim[t] = Y_val
        I_sim[t] = I_val
        ACk_sim[t] = ACk_val
        ACl_sim[t] = ACl_val
        K_sim[t + 1] = Kprime_val
        H_sim[t] = H_val
        L_sim[t] = L_val
        C_sim[t] = Y_val - I_val - ACk_val - ACl_val
    
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
        returnfirm_noise=returnfirm_noise[:, :p.nfirmspub],
        returnfirmsd_noise=returnfirmsd_noise[:, :p.nfirmspub],
        dist_zkl=dist_zkl,
        disaster_indicators=disaster_indicators,
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
    
    CRITICAL FIX: Use MT19937 matching Fortran's random_number
    
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
    
    # CRITICAL FIX: Use MT19937 matching Fortran
    rng = create_mt19937_rng(seed)
    
    # Storage for IRFs.
    irf_Y = np.zeros((T_irf, n_sims))
    irf_I = np.zeros((T_irf, n_sims))

    # Fortran controls are 1-based; convert to 0-based.
    shock_period = max(1, int(p.shockperIRF) - 1)  # 0-based
    window_start = max(1, int(p.numdiscIRF) - 1)
    # Ensure enough history + requested horizon.
    T_full = max(int(p.lengthIRF), window_start + T_irf + 1)

    for sim in range(n_sims):
        # Common random numbers across baseline and shocked paths.
        u_s = rng.random(T_full)
        u_a = rng.random(T_full)

        def _build_paths(force_shock: bool):
            a_path = np.zeros(T_full, dtype=int)
            s_path = np.zeros(T_full, dtype=int)
            a_path[0] = max(0, min(p.ainit - 1, p.anum - 1))
            s_path[0] = max(0, min(p.sinit - 1, p.snum - 1))
            for t in range(1, T_full):
                if force_shock and shock_size > 0:
                    # Match Fortran uncexogsimIRF:
                    # singleshock=1: one high-unc period at shockperIRF
                    # singleshock=0: low-unc block then one high-unc period.
                    if int(getattr(p, "singleshock", 1)) == 1:
                        if t == shock_period:
                            s_path[t] = 1
                        else:
                            cdf_s = np.cumsum(g.pr_mat_s[s_path[t - 1], :])
                            s_path[t] = int(np.searchsorted(cdf_s, u_s[t]))
                            s_path[t] = min(s_path[t], p.snum - 1)
                    else:
                        low_start = shock_period
                        low_end = shock_period + int(p.shocklengthIRF)
                        high_once = low_end + 1
                        if low_start <= t <= low_end:
                            s_path[t] = 0
                        elif t == high_once:
                            s_path[t] = 1
                        else:
                            cdf_s = np.cumsum(g.pr_mat_s[s_path[t - 1], :])
                            s_path[t] = int(np.searchsorted(cdf_s, u_s[t]))
                            s_path[t] = min(s_path[t], p.snum - 1)
                else:
                    cdf_s = np.cumsum(g.pr_mat_s[s_path[t - 1], :])
                    s_path[t] = int(np.searchsorted(cdf_s, u_s[t]))
                    s_path[t] = min(s_path[t], p.snum - 1)

                cdf_a = np.cumsum(g.pr_mat_a[a_path[t - 1], :, s_path[t - 1]])
                a_path[t] = int(np.searchsorted(cdf_a, u_a[t]))
                a_path[t] = min(a_path[t], p.anum - 1)
            return a_path, s_path

        a_base, s_base = _build_paths(force_shock=False)
        a_shock, s_shock = _build_paths(force_shock=True)

        # Same seed keeps firm-level Monte Carlo noise aligned across paths.
        sim_seed = seed + sim
        base = simulate_all_firms(
            p, g, vfi_sol, T=T_full, seed=sim_seed, verbose=False,
            a_pos_override=a_base, s_pos_override=s_base
        )
        shock = simulate_all_firms(
            p, g, vfi_sol, T=T_full, seed=sim_seed, verbose=False,
            a_pos_override=a_shock, s_pos_override=s_shock
        )

        # Figure 8 convention: percent deviation from baseline path.
        post = slice(window_start, window_start + T_irf)
        y_base = np.maximum(base.Y_sim[post], 1e-12)
        i_base_raw = base.I_sim[post]
        i_base = np.where(
            np.abs(i_base_raw) > 1e-8,
            i_base_raw,
            np.where(i_base_raw >= 0.0, 1e-8, -1e-8),
        )
        irf_Y[:, sim] = (shock.Y_sim[post] / y_base - 1.0) * 100.0
        irf_I[:, sim] = (shock.I_sim[post] / i_base - 1.0) * 100.0
    
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
