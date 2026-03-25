"""
Value Function Iteration (VFI) for BBT (2024) Model.

Implements the Bellman equation solution matching Fortran VOL_GROWTH_wrapper.f90 exactly:
1. Howard acceleration (policy iteration) - accelmaxit iterations
2. Expected value matrix computation (EVmat)
3. RHS optimization step
4. Forecast rules for aggregate capital and prices
5. Interpolation over aggregate capital forecast grid

Key arrays (matching Fortran):
- V(endog, exog, fcst): Value function
- polmat(endog, exog, fcst): Policy indices
- EVmat(exog, pol, fcst): Expected continuation values
- kbarfcstinds: Interpolation indices for aggregate capital
- kbarfcstweights: Interpolation weights
- pfcstmat: Price forecast matrix
- WLmat: Wage bill matrix
"""

import numpy as np
from typing import Tuple, Optional, Dict
from dataclasses import dataclass
from numba import jit, prange
import warnings

from .params import ModelParameters
from .grids import StateGrids, get_state_indices
from .adjustment import (
    output, capital_adjustment_cost, labor_adjustment_cost,
    AdjustmentCostCalculator
)


@dataclass
class VFISolution:
    """Container for VFI solution."""
    
    V: np.ndarray              # Value function, shape (numendog, numexog, kbarnum)
    polmat: np.ndarray         # Policy indices, shape (numendog, numexog, kbarnum)
    kprime_pos: np.ndarray     # Capital policy indices
    lpol_pos: np.ndarray       # Labor policy indices
    EVmat: np.ndarray          # Expected value matrix
    converged: bool
    iterations: int
    vf_error: float
    pol_error: float


@dataclass
class ForecastMatrices:
    """Container for forecast matrices used in VFI."""
    
    kbarfcstinds: np.ndarray    # Interpolation indices, shape (anum, snum, snum, kbarnum)
    kbarfcstweights: np.ndarray # Interpolation weights, shape (anum, snum, snum, kbarnum)
    pfcstmat: np.ndarray        # Price forecast, shape (anum, snum, snum, kbarnum)
    wfcstmat: np.ndarray        # Wage forecast, shape (anum, snum, snum, kbarnum)
    kbarfcstmat: np.ndarray     # Aggregate capital forecast, shape (anum, snum, snum, kbarnum)


def initialize_forecast_matrices(
    params: ModelParameters,
    grids: StateGrids
) -> ForecastMatrices:
    """
    Initialize forecast matrices.
    
    These are used for interpolation of aggregate capital and prices.
    
    Matches Fortran:
    ```fortran
    kbarfcstmat(:,:,:,:) = (kbar0(1)+kbar0(2))/dble(2.0)
    pfcstmat(:,:,:,:) = pval
    wfcstmat(:,:,:,:) = theta/pval
    kbarfcstinds(:,:,:,:) = 1
    kbarfcstweights(:,:,:,:) = 0.0
    ```
    """
    p = params
    g = grids
    
    # Dimensions
    anum = p.anum
    snum = p.snum
    kbarnum = p.kbarnum
    
    # Initialize with default values (Fortran lines 512-516)
    kbar_mid = (g.kbar_grid[0] + g.kbar_grid[1]) / 2.0
    
    kbarfcstmat = np.full((anum, snum, snum, kbarnum), kbar_mid)
    pfcstmat = np.full((anum, snum, snum, kbarnum), p.pval)
    wfcstmat = np.full((anum, snum, snum, kbarnum), p.theta / p.pval)
    # IMPORTANT: Fortran uses 1-based indexing, Python uses 0-based
    # Fortran initializes kbarfcstinds = 1, which corresponds to index 0 in Python
    # This is used for interpolation: weight * V[..., ind+1] + (1-weight) * V[..., ind]
    # In Fortran: ind=1, weight=0 -> V[..., 2] * 0 + V[..., 1] * 1 = V[1] (first element)
    # In Python: ind=0, weight=0 -> V[..., 1] * 0 + V[..., 0] * 1 = V[0] (first element)
    kbarfcstinds = np.zeros((anum, snum, snum, kbarnum), dtype=np.int64)  # 0 in Python = 1 in Fortran
    kbarfcstweights = np.zeros((anum, snum, snum, kbarnum))
    
    return ForecastMatrices(
        kbarfcstinds=kbarfcstinds,
        kbarfcstweights=kbarfcstweights,
        pfcstmat=pfcstmat,
        wfcstmat=wfcstmat,
        kbarfcstmat=kbarfcstmat
    )


def build_return_matrices(
    params: ModelParameters,
    grids: StateGrids,
    adj_calc: AdjustmentCostCalculator,
    fcst_mats: ForecastMatrices
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Pre-compute all return matrices for VFI.
    
    Matches Fortran lines 699-783:
    - Ymat: Output matrix
    - Imat: Investment matrix
    - ACkmat: Capital adjustment cost matrix
    - AClmat: Labor adjustment cost matrix
    - WLmat: Wage bill matrix
    
    Returns
    -------
    Ymat, Imat, ACkmat, AClmat, WLmat : np.ndarray
        Pre-computed matrices.
    """
    p = params
    g = grids
    
    znum = p.znum
    anum = p.anum
    knum = p.knum
    lnum = p.lnum
    snum = p.snum
    kbarnum = p.kbarnum
    
    # Output matrix: Ymat(z, a, k, l)
    # Fortran line 699-712
    Ymat = np.zeros((znum, anum, knum, lnum))
    for zct in range(znum):
        for act in range(anum):
            for kct in range(knum):
                for lct in range(lnum):
                    Ymat[zct, act, kct, lct] = output(
                        g.z_grid[zct], g.a_grid[act],
                        g.k_grid[kct], g.l_grid[lct],
                        p.alpha, p.nu
                    )
    
    # Investment matrix: Imat(k, k')
    # Fortran line 715-720
    Imat = np.zeros((knum, knum))
    for kct in range(knum):
        for kprimect in range(knum):
            Imat[kct, kprimect] = g.k_grid[kprimect] - (1 - p.deltak) * g.k_grid[kct]
    
    # Capital adjustment cost matrix: ACkmat(z, a, k, l, k')
    # Fortran line 723-739
    ACkmat = np.zeros((znum, anum, knum, lnum, knum))
    for zct in range(znum):
        for act in range(anum):
            for kct in range(knum):
                for lct in range(lnum):
                    for kprimect in range(knum):
                        ACkmat[zct, act, kct, lct, kprimect] = capital_adjustment_cost(
                            g.k_grid[kprimect], g.k_grid[kct],
                            p.capirrev, p.capfix, p.deltak
                        )
    
    # Labor adjustment cost matrix: AClmat(exog, k, l, l_{-1}, fcst)
    # This depends on wage which depends on forecast state
    # Fortran line 742-766
    numexog = g.numexog
    AClmat = np.zeros((numexog, knum, lnum, lnum, kbarnum))
    
    for exogct in range(numexog):
        act = g.exog_pos[exogct, 1]
        sct = g.exog_pos[exogct, 2]
        smin1ct = g.exog_pos[exogct, 3]
        
        for kct in range(knum):
            for lct in range(lnum):
                for lmin1ct in range(lnum):
                    for fcstct in range(kbarnum):
                        wfcstval = fcst_mats.wfcstmat[act, sct, smin1ct, fcstct]
                        AClmat[exogct, kct, lct, lmin1ct, fcstct] = labor_adjustment_cost(
                            g.l_grid[lct], g.l_grid[lmin1ct], wfcstval,
                            p.hirelin, p.firelin, p.labfix
                        )
    
    # Wage bill matrix: WLmat(a, s, s_{-1}, fcst, l)
    # Fortran line 769-783
    WLmat = np.zeros((anum, snum, snum, kbarnum, lnum))
    for act in range(anum):
        for sct in range(snum):
            for smin1ct in range(snum):
                for fcstct in range(kbarnum):
                    for lct in range(lnum):
                        WLmat[act, sct, smin1ct, fcstct, lct] = (
                            fcst_mats.wfcstmat[act, sct, smin1ct, fcstct] * g.l_grid[lct]
                        )
    
    return Ymat, Imat, ACkmat, AClmat, WLmat


@jit(nopython=True, parallel=True)
def howard_acceleration_step(
    V: np.ndarray,
    polmat: np.ndarray,
    Ymat: np.ndarray,
    Imat: np.ndarray,
    ACkmat: np.ndarray,
    AClmat: np.ndarray,
    WLmat: np.ndarray,
    pfcstmat: np.ndarray,
    kbarfcstinds: np.ndarray,
    kbarfcstweights: np.ndarray,
    pr_mat: np.ndarray,
    exog_pos: np.ndarray,
    endog_pos: np.ndarray,
    beta: float,
    numendog: int,
    numexog: int,
    kbarnum: int,
    numexog_next: int,
    n_accel: int,
    znum: int,
    anum: int,
    knum: int,
    lnum: int,
    snum: int
) -> np.ndarray:
    """
    Howard acceleration: evaluate Bellman with fixed policies.
    
    Matches Fortran lines 814-867 exactly.
    
    V(s) = R(s, pi(s)) + beta * E[V(s') | s, pi(s)]
    """
    V_new = V.copy()
    
    for accel in range(n_accel):
        V_old = V_new.copy()
        
        # Loop over all states
        for ct in prange(numendog * numexog * kbarnum):
            # Extract indices
            endogct = ct // (numexog * kbarnum)
            remainder = ct % (numexog * kbarnum)
            exogct = remainder // kbarnum
            fcstct = remainder % kbarnum
            
            # Extract exogenous state positions
            zct = exog_pos[exogct, 0]
            act = exog_pos[exogct, 1]
            sct = exog_pos[exogct, 2]
            smin1ct = exog_pos[exogct, 3]
            
            # Extract endogenous state positions
            kct = endog_pos[endogct, 0]
            lmin1ct = endog_pos[endogct, 1]
            
            # Get policy
            polstar = polmat[endogct, exogct, fcstct]
            
            # Extract policy positions
            kprimect = endog_pos[polstar, 0]
            lct = endog_pos[polstar, 1]
            
            # Period return (Fortran lines 837-839)
            period_return = pfcstmat[act, sct, smin1ct, fcstct] * (
                Ymat[zct, act, kct, lct]
                - ACkmat[zct, act, kct, lct, kprimect]
                - AClmat[exogct, kct, lct, lmin1ct, fcstct]
                - Imat[kct, kprimect]
                - WLmat[act, sct, smin1ct, fcstct, lct]
            )
            
            # Interpolation for aggregate capital forecast
            ind = kbarfcstinds[act, sct, smin1ct, fcstct]
            weight = kbarfcstweights[act, sct, smin1ct, fcstct]
            
            # Expected continuation value
            EV = 0.0
            for exogprimect in range(numexog_next):
                # Interpolated continuation value (Fortran lines 850-851)
                # Note: kbarnum=2, so ind can be 0 or 1 (after our fix)
                # When ind=1, ind+1=2 would be out of bounds, so clamp
                if ind >= kbarnum - 1:
                    # At upper boundary, just use the upper value
                    Vnextval = V_old[polstar, exogprimect, kbarnum - 1]
                else:
                    Vnextval = (weight * V_old[polstar, exogprimect, ind + 1] +
                               (1.0 - weight) * V_old[polstar, exogprimect, ind])
                
                EV += pr_mat[exogct, exogprimect] * Vnextval
            
            V_new[endogct, exogct, fcstct] = period_return + beta * EV
    
    return V_new


@jit(nopython=True, parallel=True)
def compute_ev_matrix(
    V: np.ndarray,
    pr_mat: np.ndarray,
    kbarfcstinds: np.ndarray,
    kbarfcstweights: np.ndarray,
    exog_pos: np.ndarray,
    numendog: int,
    numexog: int,
    kbarnum: int,
    numexog_next: int,
    anum: int,
    snum: int
) -> np.ndarray:
    """
    Compute expected continuation values for all state-policy combinations.
    
    Matches Fortran lines 870-902.
    
    EV[exog, pol, fcst] = sum over exog' of P[exog, exog'] * V[pol, exog', fcst_interpolated]
    """
    EVmat = np.zeros((numexog, numendog, kbarnum))
    
    for ct in prange(numendog * numexog * kbarnum):
        polct = ct // (numexog * kbarnum)
        remainder = ct % (numexog * kbarnum)
        exogct = remainder // kbarnum
        fcstct = remainder % kbarnum
        
        # Get act and sct
        act = exog_pos[exogct, 1]
        sct = exog_pos[exogct, 2]
        smin1ct = exog_pos[exogct, 3]
        
        # Interpolation
        ind = kbarfcstinds[act, sct, smin1ct, fcstct]
        weight = kbarfcstweights[act, sct, smin1ct, fcstct]
        
        ev = 0.0
        for exogprimect in range(numexog_next):
            # Boundary check to prevent index out of bounds
            if ind >= kbarnum - 1:
                Vnextval = V[polct, exogprimect, kbarnum - 1]
            else:
                Vnextval = (weight * V[polct, exogprimect, ind + 1] +
                           (1.0 - weight) * V[polct, exogprimect, ind])
            ev += pr_mat[exogct, exogprimect] * Vnextval
        
        EVmat[exogct, polct, fcstct] = ev
    
    return EVmat


@jit(nopython=True, parallel=True)
def optimization_step(
    Ymat: np.ndarray,
    Imat: np.ndarray,
    ACkmat: np.ndarray,
    AClmat: np.ndarray,
    WLmat: np.ndarray,
    pfcstmat: np.ndarray,
    EVmat: np.ndarray,
    exog_pos: np.ndarray,
    endog_pos: np.ndarray,
    beta: float,
    numendog: int,
    numexog: int,
    kbarnum: int,
    znum: int,
    anum: int,
    knum: int,
    lnum: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find optimal policies for each state.
    
    Matches Fortran lines 904-943.
    
    For each state, find policy that maximizes:
    RHS[pol] = R(s, pol) + beta * EV[exog, pol, fcst]
    """
    V_new = np.zeros((numendog, numexog, kbarnum))
    polmat_new = np.zeros((numendog, numexog, kbarnum), dtype=np.int64)
    
    for ct in prange(numendog * numexog * kbarnum):
        endogct = ct // (numexog * kbarnum)
        remainder = ct % (numexog * kbarnum)
        exogct = remainder // kbarnum
        fcstct = remainder % kbarnum
        
        # Extract state positions
        zct = exog_pos[exogct, 0]
        act = exog_pos[exogct, 1]
        sct = exog_pos[exogct, 2]
        smin1ct = exog_pos[exogct, 3]
        kct = endog_pos[endogct, 0]
        lmin1ct = endog_pos[endogct, 1]
        
        # Find optimal policy
        best_val = -1e30
        best_pol = 0
        
        for polct in range(numendog):
            kprimect = endog_pos[polct, 0]
            lct = endog_pos[polct, 1]
            
            # Period return
            period_return = pfcstmat[act, sct, smin1ct, fcstct] * (
                Ymat[zct, act, kct, lct]
                - ACkmat[zct, act, kct, lct, kprimect]
                - AClmat[exogct, kct, lct, lmin1ct, fcstct]
                - Imat[kct, kprimect]
                - WLmat[act, sct, smin1ct, fcstct, lct]
            )
            
            # Total value
            val = period_return + beta * EVmat[exogct, polct, fcstct]
            
            if val > best_val:
                best_val = val
                best_pol = polct
        
        V_new[endogct, exogct, fcstct] = best_val
        polmat_new[endogct, exogct, fcstct] = best_pol
    
    return V_new, polmat_new


def solve_vfi(
    params: ModelParameters,
    grids: StateGrids,
    fcst_mats: ForecastMatrices = None,
    max_iter: int = None,
    tol: float = None,
    verbose: bool = True
) -> VFISolution:
    """
    Solve firm value function via policy iteration with Howard acceleration.
    
    This implements the algorithm from VOL_GROWTH_wrapper.f90 exactly:
    
    1. Initialize V and policies
    2. For each VFI iteration:
       a. Howard acceleration: Evaluate Bellman with fixed policies (accelmaxit times)
       b. Compute EVmat for all state-policy combinations
       c. Optimization: Find optimal policies for each state
       d. Check convergence on policy function
    
    Parameters
    ----------
    params : ModelParameters
        Model parameters.
    grids : StateGrids
        State grids container.
    fcst_mats : ForecastMatrices, optional
        Forecast matrices. If None, initializes with defaults.
    max_iter : int, optional
        Max VFI iterations. Default from params.vfmaxit.
    tol : float, optional
        Convergence tolerance. Default from params.vferrortol.
    verbose : bool
        Print progress.
    
    Returns
    -------
    VFISolution
        Solution container with V, policies, and convergence info.
    """
    p = params
    g = grids
    
    # Set defaults
    if max_iter is None:
        max_iter = p.vfmaxit
    if tol is None:
        tol = p.vferrortol
    if fcst_mats is None:
        fcst_mats = initialize_forecast_matrices(params, grids)
    
    # Dimensions
    numendog = g.numendog
    numexog = g.numexog
    kbarnum = p.kbarnum
    numexog_next = g.pr_mat_full.shape[1]
    
    # Initialize (Fortran lines 795-801)
    V = np.zeros((numendog, numexog, kbarnum))
    V_old = np.zeros_like(V)
    polmat = np.zeros((numendog, numexog, kbarnum), dtype=np.int64)
    polmat_old = np.zeros_like(polmat)
    
    # Initialize policies to middle of grid
    polmat[:, :, :] = numendog // 2
    polmat_old[:, :, :] = numendog // 2
    
    # Build return matrices
    if verbose:
        print("Building return matrices...")
    
    adj_calc = AdjustmentCostCalculator(params, grids)
    Ymat, Imat, ACkmat, AClmat, WLmat = build_return_matrices(
        params, grids, adj_calc, fcst_mats
    )
    
    if verbose:
        print("Starting VFI...")
        print(f"  States: {numendog} endog x {numexog} exog x {kbarnum} fcst")
        print(f"  Max iterations: {max_iter}, Tolerance: {tol}")
        print(f"  Howard acceleration: {p.accelmaxit} steps per VFI iteration")
    
    converged = False
    vf_error = 0.0
    pol_error = 0.0
    
    for vf_iter in range(max_iter):
        
        # =====================
        # Howard Acceleration Step
        # =====================
        V = howard_acceleration_step(
            V, polmat_old,
            Ymat, Imat, ACkmat, AClmat, WLmat,
            fcst_mats.pfcstmat,
            fcst_mats.kbarfcstinds,
            fcst_mats.kbarfcstweights,
            g.pr_mat_full,
            g.exog_pos,
            g.endog_pos,
            p.beta,
            numendog, numexog, kbarnum, numexog_next, p.accelmaxit,
            p.znum, p.anum, p.knum, p.lnum, p.snum
        )
        
        # =====================
        # Compute EVmat
        # =====================
        EVmat = compute_ev_matrix(
            V, g.pr_mat_full,
            fcst_mats.kbarfcstinds,
            fcst_mats.kbarfcstweights,
            g.exog_pos,
            numendog, numexog, kbarnum, numexog_next,
            p.anum, p.snum
        )
        
        # =====================
        # Optimization Step
        # =====================
        V_new, polmat_new = optimization_step(
            Ymat, Imat, ACkmat, AClmat, WLmat,
            fcst_mats.pfcstmat,
            EVmat,
            g.exog_pos, g.endog_pos,
            p.beta,
            numendog, numexog, kbarnum,
            p.znum, p.anum, p.knum, p.lnum
        )
        
        # =====================
        # Check Convergence
        # =====================
        vf_error = np.max(np.abs(V_new - V))
        pol_error = np.max(np.abs(polmat_new - polmat_old))
        
        if verbose:
            print(f"  VFI iter {vf_iter + 1}: VF error = {vf_error:.6f}, Policy error = {pol_error:.6f}")
        
        # Update
        V = V_new
        polmat = polmat_new
        
        # Convergence check (Fortran line 957: exit on policy convergence)
        if pol_error < tol:
            converged = True
            if verbose:
                print(f"  Converged at iteration {vf_iter + 1}")
            break
        
        V_old = V.copy()
        polmat_old = polmat.copy()
    
    # Check for grid boundary hits (Fortran lines 987-1005)
    kprime_pos = np.zeros_like(polmat)
    lpol_pos = np.zeros_like(polmat)
    
    for endogct in range(numendog):
        for exogct in range(numexog):
            for fcstct in range(kbarnum):
                pol_idx = polmat[endogct, exogct, fcstct]
                kprime_pos[endogct, exogct, fcstct] = g.endog_pos[pol_idx, 0]
                lpol_pos[endogct, exogct, fcstct] = g.endog_pos[pol_idx, 1]
    
    if verbose:
        if np.max(kprime_pos) == p.knum - 1:
            print("  Warning: Hit top of capital grid")
        if np.min(kprime_pos) == 0:
            print("  Warning: Hit bottom of capital grid")
        if np.max(lpol_pos) == p.lnum - 1:
            print("  Warning: Hit top of labor grid")
        if np.min(lpol_pos) == 0:
            print("  Warning: Hit bottom of labor grid")
    
    return VFISolution(
        V=V,
        polmat=polmat,
        kprime_pos=kprime_pos,
        lpol_pos=lpol_pos,
        EVmat=EVmat,
        converged=converged,
        iterations=vf_iter + 1,
        vf_error=vf_error,
        pol_error=pol_error
    )


def solve_vfi_simplified(
    params: ModelParameters,
    grids: StateGrids,
    price: float = 1.34,
    max_iter: int = 100,
    tol: float = 1e-4,
    howard_accel: int = 20,
    verbose: bool = True
) -> VFISolution:
    """
    Simplified VFI for faster computation with Howard acceleration.
    
    This version uses fixed price and wage (no GE iteration).
    Useful for testing and IRF computation.
    
    Parameters
    ----------
    params : ModelParameters
    grids : StateGrids
    price : float
        Output price (fixed).
    max_iter : int
        Max VFI iterations (outer loop).
    tol : float
        Convergence tolerance.
    howard_accel : int
        Howard acceleration iterations per VFI iteration.
    verbose : bool
    
    Returns
    -------
    VFISolution
    """
    p = params
    g = grids
    
    # Create simple forecast matrices with fixed price
    fcst_mats = initialize_forecast_matrices(params, grids)
    fcst_mats.pfcstmat[:, :, :, :] = price
    fcst_mats.wfcstmat[:, :, :, :] = p.theta / price
    
    # Run full VFI
    return solve_vfi(params, grids, fcst_mats, max_iter, tol, verbose)
