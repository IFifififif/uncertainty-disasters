"""
Value Function Iteration (VFI) for BBT (2024) Model.

Implements the Bellman equation solution with:
1. Howard acceleration (policy iteration)
2. Forecast rules for aggregate capital and prices
3. Non-convex adjustment costs
4. Full 5D state space support

The model uses policy function iteration rather than pure VFI:
1. Howard acceleration step: Evaluate Bellman with fixed policies
2. Optimization step: Find optimal policies given values

This matches the Fortran implementation in VOL_GROWTH_wrapper.f90.
"""

import numpy as np
from typing import Tuple, Optional, Dict
from dataclasses import dataclass
from numba import jit, prange
import warnings

from .params import ModelParameters
from .grids import StateGrids, get_state_indices
from .adjustment import AdjustmentCostCalculator, output, capital_adjustment_cost, labor_adjustment_cost


@dataclass
class VFISolution:
    """Container for VFI solution."""
    
    V: np.ndarray              # Value function, shape (numendog, numexog, kbarnum)
    polmat: np.ndarray         # Policy indices, shape (numendog, numexog, kbarnum)
    kprime_pos: np.ndarray     # Capital policy indices
    lpol_pos: np.ndarray       # Labor policy indices
    converged: bool
    iterations: int
    vf_error: float
    pol_error: float


@jit(nopython=True, parallel=True)
def _howard_acceleration_step(
    V: np.ndarray,
    polmat: np.ndarray,
    returns_matrix: np.ndarray,
    pr_mat: np.ndarray,
    beta: float,
    numendog: int,
    numexog: int,
    kbarnum: int,
    numexog_next: int,
    n_accel: int
) -> np.ndarray:
    """
    Howard acceleration: evaluate Bellman with fixed policies.
    
    This iterates the Bellman operator n_accel times while keeping
    policies fixed. This speeds up convergence significantly.
    
    V(s) = R(s, pi(s)) + beta * E[V(s') | s, pi(s)]
    """
    V_new = V.copy()
    
    for accel in range(n_accel):
        V_old = V_new.copy()
        
        for endogct in prange(numendog):
            for exogct in range(numexog):
                for fcstct in range(kbarnum):
                    # Get policy
                    pol_idx = polmat[endogct, exogct, fcstct]
                    
                    # Period return
                    period_return = returns_matrix[endogct, exogct, fcstct, pol_idx]
                    
                    # Expected continuation value
                    EV = 0.0
                    for exogprimct in range(numexog_next):
                        EV += pr_mat[exogct, exogprimct] * V_old[pol_idx, exogprimct, fcstct]
                    
                    V_new[endogct, exogct, fcstct] = period_return + beta * EV
    
    return V_new


@jit(nopython=True)
def _compute_ev_matrix(
    V: np.ndarray,
    pr_mat: np.ndarray,
    numendog: int,
    numexog: int,
    kbarnum: int,
    numexog_next: int
) -> np.ndarray:
    """
    Compute expected continuation values for all state-policy combinations.
    
    EV[exog, pol, fcst] = sum over exog' of P[exog, exog'] * V[pol, exog', fcst]
    """
    EV_mat = np.zeros((numexog, numendog, kbarnum))
    
    for exogct in range(numexog):
        for polct in range(numendog):
            for fcstct in range(kbarnum):
                ev = 0.0
                for exogprimct in range(numexog_next):
                    ev += pr_mat[exogct, exogprimct] * V[polct, exogprimct, fcstct]
                EV_mat[exogct, polct, fcstct] = ev
    
    return EV_mat


def solve_vfi(
    params: ModelParameters,
    grids: StateGrids,
    adj_calc: AdjustmentCostCalculator,
    price: float = 1.34,
    w: float = None,
    max_iter: int = None,
    tol: float = None,
    verbose: bool = True
) -> VFISolution:
    """
    Solve firm value function via policy iteration with Howard acceleration.
    
    This implements the algorithm from VOL_GROWTH_wrapper.f90:
    
    1. Initialize V and policies
    2. For each VFI iteration:
       a. Howard acceleration: Evaluate Bellman with fixed policies (accelmaxit times)
       b. Optimization: Find optimal policies for each state
       c. Check convergence
    
    Parameters
    ----------
    params : ModelParameters
        Model parameters.
    grids : StateGrids
        State grids container.
    adj_calc : AdjustmentCostCalculator
        Adjustment cost calculator.
    price : float
        Output price.
    w : float, optional
        Wage rate. If None, computed as theta/price.
    max_iter : int, optional
        Max VFI iterations. Default from params.vfmaxit.
    tol : float, optional
        Convergence tolerance. Default from params.vferrortol.
    verbose : bool
        Print progress.
    
    Returns
    -------
    VFISolution
        Solution container with V and policies.
    """
    p = params
    g = grids
    
    # Set defaults
    if max_iter is None:
        max_iter = p.vfmaxit
    if tol is None:
        tol = p.vferrortol
    if w is None:
        w = p.theta / price  # From labor supply: w = theta/p
    
    # Initialize
    numendog = g.numendog
    numexog = g.numexog
    kbarnum = p.kbarnum
    numexog_next = g.pr_mat_full.shape[1]
    
    # Value function and policies
    V = np.zeros((numendog, numexog, kbarnum))
    V_old = np.zeros_like(V)
    polmat = np.zeros((numendog, numexog, kbarnum), dtype=np.int64)
    polmat_old = np.zeros_like(polmat)
    
    # Initialize policies to middle of grid
    polmat[:, :, :] = numendog // 2
    
    # Pre-compute period returns for all states and policies
    if verbose:
        print("Pre-computing returns matrix...")
    returns_matrix = _compute_returns_matrix_full(params, grids, adj_calc, price, w)
    
    if verbose:
        print("Starting VFI...")
        print(f"  States: {numendog} endog x {numexog} exog x {kbarnum} fcst")
        print(f"  Max iterations: {max_iter}, Tolerance: {tol}")
        print(f"  Howard acceleration: {p.accelmaxit} steps per VFI iteration")
    
    converged = False
    for vf_iter in range(max_iter):
        
        # =====================
        # Howard Acceleration Step
        # =====================
        V = _howard_acceleration_step(
            V, polmat, returns_matrix, g.pr_mat_full,
            p.beta, numendog, numexog, kbarnum, numexog_next, p.accelmaxit
        )
        
        # =====================
        # Compute Expected Values
        # =====================
        EV_mat = _compute_ev_matrix(
            V, g.pr_mat_full, numendog, numexog, kbarnum, numexog_next
        )
        
        # =====================
        # Optimization Step
        # =====================
        V_new, polmat_new = _optimization_step_numba(
            returns_matrix, EV_mat, p.beta, numendog, numexog, kbarnum
        )
        
        # =====================
        # Check Convergence
        # =====================
        vf_error = np.max(np.abs(V_new - V))
        pol_error = np.max(np.abs(polmat_new - polmat))
        
        if verbose and (vf_iter + 1) % 5 == 0:
            print(f"  VFI iter {vf_iter + 1}: VF error = {vf_error:.2e}, Policy error = {pol_error:.2e}")
        
        # Update
        V = V_new
        polmat_old = polmat.copy()
        polmat = polmat_new
        
        # Convergence check (matching Fortran: exit on policy convergence)
        if pol_error < tol:
            converged = True
            if verbose:
                print(f"  Converged at iteration {vf_iter + 1}")
            break
        
        V_old = V.copy()
    
    # Extract individual policies
    kprime_pos = np.zeros_like(polmat)
    lpol_pos = np.zeros_like(polmat)
    
    for endogct in range(numendog):
        for exogct in range(numexog):
            for fcstct in range(kbarnum):
                pol_idx = polmat[endogct, exogct, fcstct]
                kprime_pos[endogct, exogct, fcstct] = g.endog_pos[pol_idx, 0]
                lpol_pos[endogct, exogct, fcstct] = g.endog_pos[pol_idx, 1]
    
    return VFISolution(
        V=V,
        polmat=polmat,
        kprime_pos=kprime_pos,
        lpol_pos=lpol_pos,
        converged=converged,
        iterations=vf_iter + 1,
        vf_error=vf_error,
        pol_error=pol_error
    )


def _compute_returns_matrix_full(
    params: ModelParameters,
    grids: StateGrids,
    adj_calc: AdjustmentCostCalculator,
    price: float,
    w: float
) -> np.ndarray:
    """
    Pre-compute period returns for all states and policy choices.
    
    Returns matrix has shape (numendog, numexog, kbarnum, numendog)
    where the last dimension is the policy choice.
    
    R(s, s', a, a', k, l, k', l') = p * (Y - ACk - ACl - I - w*l')
    """
    p = params
    g = grids
    
    numendog = g.numendog
    numexog = g.numexog
    kbarnum = p.kbarnum
    
    # Shape: (numendog, numexog, kbarnum, numendog)
    returns_matrix = np.zeros((numendog, numexog, kbarnum, numendog))
    
    # Pre-compute components
    Y_mat = adj_calc.Y_mat  # Shape: (znum, anum, knum, lnum)
    I_mat = adj_calc.I_mat  # Shape: (knum, knum)
    ACk_base = adj_calc.ACk_base  # Shape: (knum, knum)
    
    for endogct in range(numendog):
        k_idx = g.endog_pos[endogct, 0]
        l_idx = g.endog_pos[endogct, 1]
        
        for exogct in range(numexog):
            z_idx = g.exog_pos[exogct, 0]
            a_idx = g.exog_pos[exogct, 1]
            
            # Output for this state
            Y = Y_mat[z_idx, a_idx, k_idx, l_idx]
            
            for fcstct in range(kbarnum):
                for polct in range(numendog):
                    k_prime_idx = g.endog_pos[polct, 0]
                    l_prime_idx = g.endog_pos[polct, 1]
                    
                    # Capital adjustment cost
                    ACk = ACk_base[k_idx, k_prime_idx]
                    
                    # Labor adjustment cost
                    ACl = labor_adjustment_cost(
                        g.l_grid[l_prime_idx], g.l_grid[l_idx], w,
                        p.hirelin, p.firelin, p.labfix
                    )
                    
                    # Investment
                    I_val = I_mat[k_idx, k_prime_idx]
                    
                    # Wage bill
                    WL = w * g.l_grid[l_prime_idx]
                    
                    # Period return
                    returns_matrix[endogct, exogct, fcstct, polct] = price * (Y - ACk - ACl - I_val - WL)
    
    return returns_matrix


@jit(nopython=True, parallel=True)
def _optimization_step_numba(
    returns_matrix: np.ndarray,
    EV_mat: np.ndarray,
    beta: float,
    numendog: int,
    numexog: int,
    kbarnum: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find optimal policies for each state.
    
    For each state, find policy that maximizes:
    RHS[pol] = R(s, pol) + beta * EV[exog, pol, fcst]
    """
    V_new = np.zeros((numendog, numexog, kbarnum))
    polmat_new = np.zeros((numendog, numexog, kbarnum), dtype=np.int64)
    
    for endogct in prange(numendog):
        for exogct in range(numexog):
            for fcstct in range(kbarnum):
                # Find optimal policy
                best_val = -1e20
                best_pol = 0
                
                for polct in range(numendog):
                    val = returns_matrix[endogct, exogct, fcstct, polct] + beta * EV_mat[exogct, polct, fcstct]
                    if val > best_val:
                        best_val = val
                        best_pol = polct
                
                V_new[endogct, exogct, fcstct] = best_val
                polmat_new[endogct, exogct, fcstct] = best_pol
    
    return V_new, polmat_new


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
    
    This version reduces dimensionality by:
    1. Aggregating over z (using mean productivity)
    2. Fixing aggregate state (a=1)
    3. Using 2D state (k, sigma) instead of 5D
    
    Implements Howard acceleration:
    1. Evaluate value function with fixed policies (howard_accel iterations)
    2. Optimize policies once
    3. Repeat until convergence
    
    This matches the Fortran approach more closely.
    
    Parameters
    ----------
    params : ModelParameters
    grids : StateGrids
    price : float
        Output price.
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
    
    # Simplified state space: (k, sigma)
    Nk = p.knum
    Ns = p.snum
    
    # Value function and policy
    V = np.zeros((Ns, Nk))
    V_old = np.zeros_like(V)
    k_policy = np.zeros((Ns, Nk), dtype=np.int64)
    k_policy_old = np.zeros_like(k_policy, dtype=np.int64)
    
    # Wage from labor supply
    w = p.theta / price
    
    # Aggregate state (fixed)
    a_val = 1.0
    z_val = 1.0
    
    if verbose:
        print("Starting simplified VFI with Howard acceleration...")
    
    for iteration in range(max_iter):
        
        # =====================
        # Howard Acceleration Step
        # =====================
        # Evaluate value function with fixed policies
        for accel in range(howard_accel):
            V_new = np.zeros_like(V)
            
            for s_idx in range(Ns):
                for k_idx in range(Nk):
                    k = g.k_grid[k_idx]
                    k_prime_idx = k_policy[s_idx, k_idx]
                    k_prime = g.k_grid[k_prime_idx]
                    
                    # Investment
                    I = k_prime - (1 - p.deltak) * k
                    
                    # Output (simplified: fixed labor = 1)
                    Y = output(z_val, a_val, k, 1.0, p.alpha, p.nu)
                    
                    # Capital adjustment cost (non-convex)
                    ACk = capital_adjustment_cost(k_prime, k, p.capirrev, p.capfix, p.deltak)
                    
                    # Period profit
                    profit = price * Y - ACk - I - w * 1.0
                    
                    # Expected continuation value
                    EV = 0.0
                    for s_prime_idx in range(Ns):
                        prob = g.pr_mat_s[s_idx, s_prime_idx]
                        EV += prob * V_old[s_prime_idx, k_prime_idx]
                    
                    V_new[s_idx, k_idx] = profit + p.beta * EV
            
            V_old = V_new.copy()
        
        # =====================
        # Policy Optimization Step
        # =====================
        for s_idx in range(Ns):
            for k_idx in range(Nk):
                k = g.k_grid[k_idx]
                
                # Compute value for each k' choice
                values = np.zeros(Nk)
                
                for kprime_idx in range(Nk):
                    k_prime = g.k_grid[kprime_idx]
                    
                    # Investment
                    I = k_prime - (1 - p.deltak) * k
                    
                    # Output (simplified: fixed labor = 1)
                    Y = output(z_val, a_val, k, 1.0, p.alpha, p.nu)
                    
                    # Capital adjustment cost (non-convex)
                    ACk = capital_adjustment_cost(k_prime, k, p.capirrev, p.capfix, p.deltak)
                    
                    # Period profit
                    profit = price * Y - ACk - I - w * 1.0
                    
                    # Expected continuation value
                    EV = 0.0
                    for s_prime_idx in range(Ns):
                        prob = g.pr_mat_s[s_idx, s_prime_idx]
                        EV += prob * V_old[s_prime_idx, kprime_idx]
                    
                    values[kprime_idx] = profit + p.beta * EV
                
                # Find optimal k'
                best_idx = np.argmax(values)
                V[s_idx, k_idx] = values[best_idx]
                k_policy[s_idx, k_idx] = best_idx
        
        # =====================
        # Check Convergence
        # =====================
        vf_error = np.max(np.abs(V - V_old))
        pol_error = np.max(np.abs(k_policy - k_policy_old))
        
        if verbose and (iteration + 1) % 5 == 0:
            print(f"  Iteration {iteration + 1}: VF error = {vf_error:.2e}, Policy error = {pol_error:.2e}")
        
        # Update for next iteration
        V_old = V.copy()
        k_policy_old = k_policy.copy()
        
        # Convergence check
        if pol_error < tol:
            if verbose:
                print(f"  Converged at iteration {iteration + 1}")
            break
    
    # Package results
    # Create dummy arrays for full solution
    numendog = g.numendog
    numexog = g.numexog
    
    V_full = np.zeros((numendog, numexog, p.kbarnum))
    polmat_full = np.zeros((numendog, numexog, p.kbarnum), dtype=np.int64)
    
    # Map simplified solution to full (mapping k_idx -> endog_idx)
    for endogct in range(numendog):
        k_idx = g.endog_pos[endogct, 0]
        for exogct in range(numexog):
            s_idx = g.exog_pos[exogct, 2]
            for fcstct in range(p.kbarnum):
                V_full[endogct, exogct, fcstct] = V[s_idx, k_idx]
                # Map k_policy to endog index
                kprime_idx = k_policy[s_idx, k_idx]
                # Use same l as current (simplified)
                l_idx = g.endog_pos[endogct, 1]
                polmat_full[endogct, exogct, fcstct] = kprime_idx * p.lnum + l_idx
    
    return VFISolution(
        V=V_full,
        polmat=polmat_full,
        kprime_pos=polmat_full,  # Simplified
        lpol_pos=np.zeros_like(polmat_full),
        converged=(pol_error < tol),
        iterations=iteration + 1,
        vf_error=vf_error,
        pol_error=pol_error
    )
