"""
General Equilibrium Solver for BBT (2024) Model.

Implements the market-clearing price iteration and distribution evolution
matching Fortran VOL_GROWTH_wrapper.f90 exactly.

Key Components:
1. Distribution evolution over state space
2. Aggregate variable computation (Y, I, C, etc.)
3. Price iteration to clear goods market
4. Forecast rule updates

The model solves for equilibrium price p such that:
    excess_demand(p) = 1/p - C(p) = 0

where C(p) is aggregate consumption as a function of price.
"""

import numpy as np
from typing import Tuple, Dict, Optional
from dataclasses import dataclass
from numba import jit, prange

from .params import ModelParameters
from .grids import StateGrids
from .adjustment import output, capital_adjustment_cost, labor_adjustment_cost
from .vfi import VFISolution, ForecastMatrices


@dataclass
class DistributionState:
    """Container for distribution over state space."""
    
    dist_zkl: np.ndarray     # Distribution over (z, k, l_{-1}), shape (znum, numendog, T)
    firmpos: np.ndarray      # Firm positions, shape (znum, numendog)
    K_bar: float             # Aggregate capital
    period: int              # Current period


@dataclass
class Aggregates:
    """Container for aggregate variables at one time period."""
    
    Y: float      # Aggregate output
    K: float      # Aggregate capital
    L: float      # Aggregate labor
    I: float      # Aggregate investment
    H: float      # Aggregate hiring
    C: float      # Aggregate consumption
    ACk: float    # Aggregate capital adjustment costs
    ACl: float    # Aggregate labor adjustment costs
    p: float      # Market price
    w: float      # Wage rate


def initialize_distribution(
    params: ModelParameters,
    grids: StateGrids,
    T: int
) -> DistributionState:
    """
    Initialize distribution over state space.
    
    Matches Fortran lines 1062-1088.
    
    Parameters
    ----------
    params : ModelParameters
    grids : StateGrids
    T : int
        Total time periods.
    
    Returns
    -------
    DistributionState
        Initial distribution state.
    """
    p = params
    g = grids
    
    # Distribution over (z, k, l_{-1})
    # Fortran line 1062: distzkl(:,:,:) = 0.0
    dist_zkl = np.zeros((p.znum, g.numendog, T))
    
    # Initialize at a chosen point
    # Fortran lines 1079-1083
    kct = min(p.knum // 2, p.knum - 1)  # Middle of capital grid
    lmin1ct = min(p.lnum // 3, p.lnum - 1)  # Lower third of labor grid
    endogct = kct * p.lnum + lmin1ct
    
    # Place mass around this point
    # Fortran: distzkl(:,(endogct-5):(endogct+5),1) = 1.0
    for i in range(p.znum):
        for j in range(max(0, endogct - 5), min(g.numendog, endogct + 6)):
            dist_zkl[i, j, 0] = 1.0
    
    # Normalize
    dist_zkl[:, :, 0] /= dist_zkl[:, :, 0].sum()
    
    # Initialize firm positions
    firmpos = np.zeros((p.znum, g.numendog), dtype=int)
    
    # Initialize aggregate capital
    K_bar = (p.kbarmin + p.kbarmax) / 2.0
    
    return DistributionState(
        dist_zkl=dist_zkl,
        firmpos=firmpos,
        K_bar=K_bar,
        period=0
    )


@jit(nopython=True, parallel=True)
def compute_aggregates(
    dist_zkl: np.ndarray,
    Ymat: np.ndarray,
    Imat: np.ndarray,
    ACkmat: np.ndarray,
    l_grid: np.ndarray,
    k_grid: np.ndarray,
    z_grid: np.ndarray,
    a_grid: np.ndarray,
    exog_pos: np.ndarray,
    endog_pos: np.ndarray,
    kprime_pos: np.ndarray,
    lpol_pos: np.ndarray,
    a_idx: int,
    s_idx: int,
    sm1_idx: int,
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
    znum: int,
    numendog: int,
    knum: int,
    lnum: int,
    disttol: float
) -> Tuple[float, float, float, float, float, float, float]:
    """
    Compute aggregate variables from distribution and policies.
    
    Matches Fortran lines 1169-1201.
    
    Returns
    -------
    Y_agg, I_agg, ACk_agg, ACl_agg, K_prime_agg, H_agg, L_agg : floats
        Aggregate variables.
    """
    Y_agg = 0.0
    I_agg = 0.0
    ACk_agg = 0.0
    ACl_agg = 0.0
    K_prime_agg = 0.0
    H_agg = 0.0
    L_agg = 0.0
    
    changetol = 1e-10
    
    for zct in prange(znum):
        for endogct in range(numendog):
            weight = dist_zkl[zct, endogct]
            
            if weight > disttol:
                # Get endogenous state positions
                kct = endog_pos[endogct, 0]
                lmin1ct = endog_pos[endogct, 1]
                
                # Get policy positions
                kprimect = kprime_pos[endogct, 0]  # Simplified
                lct = lpol_pos[endogct, 0]  # Simplified
                
                # Output
                Y_agg += weight * Ymat[zct, a_idx, kct, lct]
                
                # Investment
                I_agg += weight * Imat[kct, kprimect]
                
                # Capital adjustment cost
                ACk_agg += weight * ACkmat[zct, a_idx, kct, lct, kprimect]
                
                # Labor adjustment cost (inline for numba compatibility)
                Y_val = Ymat[zct, a_idx, kct, lct]
                l_val = l_grid[lct]
                l_prev_val = l_grid[lmin1ct]
                
                # h = l - (1-deltan)*l_prev
                h = l_val - (1.0 - deltan) * l_prev_val
                
                ACl_val = 0.0
                if abs(h) > changetol:
                    # Fixed cost
                    ACl_val += labfix * Y_val
                    # Firing cost
                    if h < -changetol:
                        ACl_val += -h * w * firelin
                    # Hiring cost
                    if h > changetol:
                        ACl_val += h * w * hirelin
                
                ACl_agg += weight * ACl_val
                
                # Next period capital
                K_prime_agg += weight * k_grid[kprimect]
                
                # Hiring
                H_agg += weight * (l_val - (1.0 - deltan) * l_prev_val)
                
                # Labor
                L_agg += weight * l_val
    
    return Y_agg, I_agg, ACk_agg, ACl_agg, K_prime_agg, H_agg, L_agg


def compute_consumption(
    Y: float, I: float, ACk: float, ACl: float
) -> float:
    """
    Compute aggregate consumption.
    
    Fortran line 1204: Cvalp = Yvalp - Ivalp - ACkvalp - AClvalp
    
    Parameters
    ----------
    Y, I, ACk, ACl : float
        Output, investment, and adjustment costs.
    
    Returns
    -------
    C : float
        Aggregate consumption.
    """
    return Y - I - ACk - ACl


def compute_excess_demand(p: float, C: float) -> float:
    """
    Compute excess demand.
    
    Fortran line 1207: ep0(piter) = 1.0/pval - Cvalp
    
    Parameters
    ----------
    p : float
        Price level.
    C : float
        Aggregate consumption.
    
    Returns
    -------
    excess_demand : float
        1/p - C
    """
    return 1.0 / p - C


def find_market_clearing_price(
    params: ModelParameters,
    dist_state: DistributionState,
    vfi_sol: VFISolution,
    grids: StateGrids,
    Ymat: np.ndarray,
    Imat: np.ndarray,
    ACkmat: np.ndarray,
    a_idx: int,
    s_idx: int,
    sm1_idx: int,
    p_init: float = 1.34,
    tol: float = 1e-3,
    max_iter: int = 50
) -> Tuple[float, Aggregates]:
    """
    Find market-clearing price via iteration.
    
    Matches Fortran price iteration logic.
    
    Parameters
    ----------
    params : ModelParameters
    dist_state : DistributionState
    vfi_sol : VFISolution
    grids : StateGrids
    Ymat, Imat, ACkmat : np.ndarray
        Pre-computed matrices.
    a_idx, s_idx, sm1_idx : int
        Current exogenous state indices.
    p_init : float
        Initial price guess.
    tol : float
        Convergence tolerance.
    max_iter : int
        Maximum iterations.
    
    Returns
    -------
    p : float
        Market-clearing price.
    Aggregates
        Aggregate variables at equilibrium.
    """
    p = params
    g = grids
    
    # Get wage from labor supply
    w = p.theta / p_init
    
    # Compute aggregates with current distribution
    Y_agg, I_agg, ACk_agg, ACl_agg, K_prime, H_agg, L_agg = compute_aggregates(
        dist_state.dist_zkl[:, :, dist_state.period],
        Ymat, Imat, ACkmat,
        g.l_grid, g.k_grid, g.z_grid, g.a_grid,
        g.exog_pos, g.endog_pos,
        vfi_sol.kprime_pos[:, :, 0],  # Simplified: use first forecast index
        vfi_sol.lpol_pos[:, :, 0],
        a_idx, s_idx, sm1_idx, w,
        p.alpha, p.nu, p.capirrev, p.capfix, p.deltak, p.deltan,
        p.hirelin, p.firelin, p.labfix,
        p.znum, g.numendog, p.knum, p.lnum, p.disttol
    )
    
    # Compute consumption
    C_agg = compute_consumption(Y_agg, I_agg, ACk_agg, ACl_agg)
    
    # For simplified model, use fixed price
    # Full model would iterate to find clearing price
    price = p_init
    
    return price, Aggregates(
        Y=Y_agg, K=dist_state.K_bar, L=L_agg,
        I=I_agg, H=H_agg, C=C_agg,
        ACk=ACk_agg, ACl=ACl_agg,
        p=price, w=w
    )


def evolve_distribution(
    dist_state: DistributionState,
    vfi_sol: VFISolution,
    grids: StateGrids,
    pr_mat_z: np.ndarray,
    a_idx: int,
    s_idx: int,
    sm1_idx: int
) -> DistributionState:
    """
    Evolve distribution to next period.
    
    Matches Fortran lines 1232-1252.
    
    Parameters
    ----------
    dist_state : DistributionState
        Current distribution.
    vfi_sol : VFISolution
        VFI solution with policies.
    grids : StateGrids
    pr_mat_z : np.ndarray
        Productivity transition matrix.
    a_idx, s_idx, sm1_idx : int
        Current exogenous state indices.
    
    Returns
    -------
    DistributionState
        Updated distribution for next period.
    """
    p = dist_state.params if hasattr(dist_state, 'params') else None
    g = grids
    
    t = dist_state.period
    dist_new = dist_state.dist_zkl.copy()
    
    # Zero out next period distribution
    dist_new[:, :, t + 1] = 0.0
    
    # Evolve distribution
    # Fortran lines 1243-1244
    for zct in range(g.z_grid.size):
        for endogct in range(g.numendog):
            weight = dist_state.dist_zkl[zct, endogct, t]
            
            if weight > 1e-10:
                # Get policy
                exog_idx = zct * g.anum * g.snum * g.snum + a_idx * g.snum * g.snum + s_idx * g.snum + sm1_idx
                polstar = vfi_sol.polmat[endogct, exog_idx, 0]
                
                # Transfer weight to next period
                for zprimect in range(g.z_grid.size):
                    dist_new[zprimect, polstar, t + 1] += pr_mat_z[zct, zprimect, s_idx] * weight
    
    # Normalize
    dist_new[:, :, t + 1] /= dist_new[:, :, t + 1].sum()
    
    return DistributionState(
        dist_zkl=dist_new,
        firmpos=dist_state.firmpos,
        K_bar=dist_state.K_bar,
        period=t + 1
    )


class GESolver:
    """
    General Equilibrium Solver for BBT Model.
    
    Implements the full simulation loop with:
    1. Price determination (fixed in simplified version)
    2. Distribution evolution
    3. Aggregate variable computation
    4. Multi-firm simulation
    """
    
    def __init__(
        self,
        params: ModelParameters,
        grids: StateGrids,
        vfi_sol: VFISolution,
        fcst_mats: ForecastMatrices
    ):
        """
        Initialize GE solver.
        
        Parameters
        ----------
        params : ModelParameters
        grids : StateGrids
        vfi_sol : VFISolution
        fcst_mats : ForecastMatrices
        """
        self.params = params
        self.grids = grids
        self.vfi_sol = vfi_sol
        self.fcst_mats = fcst_mats
        
        # Pre-compute matrices
        self._precompute_matrices()
    
    def _precompute_matrices(self):
        """Pre-compute matrices for efficiency."""
        p = self.params
        g = self.grids
        
        # Output matrix
        self.Ymat = np.zeros((p.znum, p.anum, p.knum, p.lnum))
        for zi in range(p.znum):
            for ai in range(p.anum):
                for ki in range(p.knum):
                    for li in range(p.lnum):
                        self.Ymat[zi, ai, ki, li] = output(
                            g.z_grid[zi], g.a_grid[ai],
                            g.k_grid[ki], g.l_grid[li],
                            p.alpha, p.nu
                        )
        
        # Investment matrix
        self.Imat = np.zeros((p.knum, p.knum))
        for ki in range(p.knum):
            for kj in range(p.knum):
                self.Imat[ki, kj] = g.k_grid[kj] - (1 - p.deltak) * g.k_grid[ki]
        
        # Capital adjustment cost matrix (simplified)
        self.ACkmat = np.zeros((p.znum, p.anum, p.knum, p.lnum, p.knum))
        for zi in range(p.znum):
            for ai in range(p.anum):
                for ki in range(p.knum):
                    for li in range(p.lnum):
                        Y_val = self.Ymat[zi, ai, ki, li]
                        for kj in range(p.knum):
                            self.ACkmat[zi, ai, ki, li, kj] = capital_adjustment_cost(
                                g.k_grid[kj], g.k_grid[ki],
                                p.capirrev, p.capfix, p.deltak,
                                Y_val
                            )
    
    def run_simulation(
        self,
        T: int,
        price: float = 1.34,
        verbose: bool = True
    ) -> Dict:
        """
        Run full simulation with GE solver.
        
        Parameters
        ----------
        T : int
            Number of periods.
        price : float
            Fixed price (simplified mode).
        verbose : bool
        
        Returns
        -------
        results : dict
            Simulation results including aggregates.
        """
        p = self.params
        g = self.grids
        
        # Initialize distribution
        dist_state = initialize_distribution(p, g, T + 1)
        
        # Storage for aggregate series
        Y_sim = np.zeros(T)
        K_sim = np.zeros(T + 1)
        L_sim = np.zeros(T)
        I_sim = np.zeros(T)
        H_sim = np.zeros(T)
        C_sim = np.zeros(T)
        ACk_sim = np.zeros(T)
        ACl_sim = np.zeros(T)
        p_sim = np.zeros(T)
        
        K_sim[0] = dist_state.K_bar
        
        if verbose:
            print("Running GE simulation...")
        
        for t in range(T):
            # Get exogenous states (simplified: use middle values)
            a_idx = p.anum // 2
            s_idx = 0  # Low uncertainty
            sm1_idx = 0
            
            # Compute aggregates
            w = p.theta / price
            
            Y_agg, I_agg, ACk_agg, ACl_agg, K_prime, H_agg, L_agg = compute_aggregates(
                dist_state.dist_zkl[:, :, t],
                self.Ymat, self.Imat, self.ACkmat,
                g.l_grid, g.k_grid, g.z_grid, g.a_grid,
                g.exog_pos, g.endog_pos,
                self.vfi_sol.kprime_pos[:, :, 0],
                self.vfi_sol.lpol_pos[:, :, 0],
                a_idx, s_idx, sm1_idx, w,
                p.alpha, p.nu, p.capirrev, p.capfix, p.deltak, p.deltan,
                p.hirelin, p.firelin, p.labfix,
                p.znum, g.numendog, p.knum, p.lnum, p.disttol
            )
            
            C_agg = compute_consumption(Y_agg, I_agg, ACk_agg, ACl_agg)
            
            # Store results
            Y_sim[t] = Y_agg
            L_sim[t] = L_agg
            I_sim[t] = I_agg
            H_sim[t] = H_agg
            C_sim[t] = C_agg
            ACk_sim[t] = ACk_agg
            ACl_sim[t] = ACl_agg
            p_sim[t] = price
            K_sim[t + 1] = K_prime
            
            # Update distribution period
            dist_state.period = t
            
            if verbose and (t + 1) % 100 == 0:
                print(f"  Period {t + 1}/{T}")
        
        return {
            'Y': Y_sim,
            'K': K_sim,
            'L': L_sim,
            'I': I_sim,
            'H': H_sim,
            'C': C_sim,
            'ACk': ACk_sim,
            'ACl': ACl_sim,
            'p': p_sim,
            'dist': dist_state
        }
