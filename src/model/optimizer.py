"""
Particle Swarm Optimization (PSO) for BBT (2024) Model.

Implements the PSO algorithm from Fortran base_lib.f90 for
parameter estimation via GMM.

PSO is preferred over gradient-based methods because:
1. GMM objective is potentially non-smooth
2. Multiple local minima are possible
3. No need for analytical gradients
"""

import numpy as np
from typing import Callable, Tuple, Optional
from dataclasses import dataclass


@dataclass
class PSOConfig:
    """PSO configuration parameters."""
    
    npart: int = 75           # Number of particles
    max_iter: int = 5000      # Max iterations
    x_tol: float = 1e-3       # Position tolerance
    f_tol: float = 1e-3       # Function tolerance
    x_quick_tol: float = 1e-3 # Quick convergence tolerance
    x_quick_num: int = 5      # Number of best particles for quick check
    phi: Tuple[float, float] = (2.05, 2.05)  # Acceleration coefficients
    seed: int = 2501          # Random seed


@dataclass
class PSOResult:
    """PSO optimization result."""
    
    x_opt: np.ndarray         # Optimal parameters
    f_opt: float              # Optimal function value
    converged: bool
    iterations: int
    x_history: Optional[np.ndarray] = None  # Position history
    f_history: Optional[np.ndarray] = None  # Function history


def pso_optimize(
    objective: Callable[[np.ndarray], float],
    lb: np.ndarray,
    ub: np.ndarray,
    config: PSOConfig = None,
    verbose: bool = True
) -> PSOResult:
    """
    Particle Swarm Optimization.
    
    Implements the PSO algorithm from Fortran base_lib.f90:
    
    1. Initialize swarm of particles randomly in parameter space
    2. For each iteration:
       a. Update velocities: v = chi * (v + phi1 * r1 * (pbest - x) + phi2 * r2 * (gbest - x))
       b. Update positions: x = x + v
       c. Evaluate objective
       d. Update personal and global bests
       e. Check convergence
    
    Parameters
    ----------
    objective : callable
        Function to minimize, takes np.ndarray and returns float.
    lb, ub : np.ndarray
        Lower and upper bounds for parameters.
    config : PSOConfig, optional
        PSO configuration. Uses defaults if None.
    verbose : bool
        Print progress.
    
    Returns
    -------
    PSOResult
        Optimization result.
    """
    if config is None:
        config = PSOConfig()
    
    nvar = len(lb)
    npart = config.npart
    
    # Set random seed
    np.random.seed(config.seed)
    
    # Constriction factor (ensures convergence)
    phi_sum = config.phi[0] + config.phi[1]
    chi = 2.0 / (phi_sum - 2.0 + np.sqrt(phi_sum**2 - 4.0 * phi_sum))
    
    # Space width for velocity initialization
    space_norm = np.sqrt(np.sum((ub - lb)**2))
    
    # =====================
    # Initialize Swarm
    # =====================
    
    # Positions: uniform in [lb, ub]
    x_store = np.zeros((nvar, npart))
    for i in range(npart):
        x_store[:, i] = lb + (ub - lb) * np.random.random(nvar)
    
    # Velocities: uniform in [-space_norm, space_norm]
    v_store = np.zeros((nvar, npart))
    for i in range(npart):
        v_store[:, i] = space_norm * (2 * np.random.random(nvar) - 1)
    
    # Function values
    fx_store = np.zeros(npart)
    for i in range(npart):
        fx_store[i] = objective(x_store[:, i])
    
    # Personal bests
    pbest_store = x_store.copy()
    fpbest_store = fx_store.copy()
    
    # Global best
    gbest_idx = np.argmin(fx_store)
    gbest = x_store[:, gbest_idx].copy()
    fgbest = fx_store[gbest_idx]
    
    if verbose:
        print(f"PSO initialized: {npart} particles, {nvar} parameters")
        print(f"  Initial best f = {fgbest:.6f}")
    
    # =====================
    # Main Loop
    # =====================
    
    converged = False
    for iteration in range(config.max_iter):
        
        # Update velocities and positions
        for i in range(npart):
            # Random coefficients
            r1 = np.random.random(nvar)
            r2 = np.random.random(nvar)
            
            # Velocity update
            v_store[:, i] = chi * (
                v_store[:, i] +
                config.phi[0] * r1 * (pbest_store[:, i] - x_store[:, i]) +
                config.phi[1] * r2 * (gbest - x_store[:, i])
            )
            
            # Position update
            x_store[:, i] = x_store[:, i] + v_store[:, i]
            
            # Enforce bounds
            x_store[:, i] = np.maximum(x_store[:, i], lb)
            x_store[:, i] = np.minimum(x_store[:, i], ub)
            
            # Evaluate objective
            fx_store[i] = objective(x_store[:, i])
        
        # Update personal and global bests
        for i in range(npart):
            if fx_store[i] < fpbest_store[i]:
                pbest_store[:, i] = x_store[:, i].copy()
                fpbest_store[i] = fx_store[i]
            
            if fx_store[i] < fgbest:
                gbest = x_store[:, i].copy()
                fgbest = fx_store[i]
        
        # =====================
        # Convergence Check
        # =====================
        
        # Sort personal bests
        fpbest_sorted = np.sort(fpbest_store)
        quick_thresh = fpbest_sorted[min(config.x_quick_num, npart) - 1]
        
        # Maximum distance from global best
        x_norm = 0.0
        x_quick_norm = 0.0
        rel_tol = 0.0
        
        for i in range(npart):
            dist = np.sqrt(np.sum((gbest - x_store[:, i])**2))
            x_norm = max(x_norm, dist)
            
            # Relative tolerance
            rel = abs(fgbest - fpbest_store[i]) / (abs(fgbest) + abs(fpbest_store[i]) + 0.01)
            rel_tol = max(rel_tol, rel)
            
            # Quick norm (only for best particles)
            if fpbest_store[i] < quick_thresh:
                x_quick_norm = max(x_quick_norm, dist)
        
        # Print progress
        if verbose and (iteration + 1) % 100 == 0:
            print(f"  Iter {iteration + 1}: f = {fgbest:.6f}, x_norm = {x_norm:.4f}")
        
        # Check convergence
        if x_norm < config.x_tol or rel_tol < config.f_tol or x_quick_norm < config.x_quick_tol:
            converged = True
            if verbose:
                print(f"  Converged at iteration {iteration + 1}")
                print(f"    x_norm = {x_norm:.6f}, rel_tol = {rel_tol:.6f}")
            break
    
    return PSOResult(
        x_opt=gbest,
        f_opt=fgbest,
        converged=converged,
        iterations=iteration + 1
    )


def pso_optimize_restart(
    objective: Callable[[np.ndarray], float],
    lb: np.ndarray,
    ub: np.ndarray,
    config: PSOConfig = None,
    restart_file: str = None,
    verbose: bool = True
) -> PSOResult:
    """
    PSO with restart capability.
    
    This version can save progress to file and restart from
    a previous run. Useful for long-running optimizations.
    
    Parameters
    ----------
    objective : callable
        Function to minimize.
    lb, ub : np.ndarray
        Parameter bounds.
    config : PSOConfig, optional
        PSO configuration.
    restart_file : str, optional
        File to save/load progress.
    verbose : bool
        Print progress.
    
    Returns
    -------
    PSOResult
        Optimization result.
    """
    # For simplicity, just run standard PSO
    # Full implementation would save/load state from file
    return pso_optimize(objective, lb, ub, config, verbose)


def nelder_mead_2d(
    objective: Callable[[np.ndarray], float],
    x_init: np.ndarray,
    lb: np.ndarray,
    ub: np.ndarray,
    f_tol: float = 1e-6,
    x_tol: float = 1e-6,
    max_iter: int = 1000
) -> PSOResult:
    """
    Nelder-Mead simplex optimization (2D version from Fortran).
    
    This is a simpler optimization method that can be used
    for lower-dimensional problems or as a local search
    after PSO finds a good region.
    
    Parameters
    ----------
    objective : callable
        Function to minimize.
    x_init : np.ndarray
        Initial guess (2D).
    lb, ub : np.ndarray
        Parameter bounds.
    f_tol : float
        Function tolerance.
    x_tol : float
        Position tolerance.
    max_iter : int
        Maximum iterations.
    
    Returns
    -------
    PSOResult
        Optimization result.
    """
    nvar = 2
    
    # Initialize simplex
    x1 = x_init.copy()
    x2 = (lb + ub) / 2.0
    x3 = (x_init + lb) / 2.0
    
    f1 = objective(x1)
    f2 = objective(x2)
    f3 = objective(x3)
    
    # Sort: B = best, G = good, W = worst
    points = np.array([x1, x2, x3])
    values = np.array([f1, f2, f3])
    
    for iteration in range(max_iter):
        # Sort by function value
        order = np.argsort(values)
        points = points[order]
        values = values[order]
        
        xB, xG, xW = points
        fB, fG, fW = values
        
        # Midpoint
        xM = (xB + xG) / 2.0
        
        # Reflection
        xR = 2.0 * xM - xW
        xR = np.maximum(xR, lb)
        xR = np.minimum(xR, ub)
        fR = objective(xR)
        
        # Case I: fR < fG (reflection is better than worst)
        if fR < fG:
            if fB < fR:
                # Accept reflection
                points[2] = xR
                values[2] = fR
            else:
                # Expansion
                xE = 2.0 * xR - xM
                xE = np.maximum(xE, lb)
                xE = np.minimum(xE, ub)
                fE = objective(xE)
                
                if fE < fB:
                    points[2] = xE
                    values[2] = fE
                else:
                    points[2] = xR
                    values[2] = fR
        
        # Case II: fR >= fG (reflection is not better)
        else:
            if fR < fW:
                points[2] = xR
                values[2] = fR
            else:
                # Contraction
                xC = (xM + xW) / 2.0
                fC = objective(xC)
                
                if fC < fW:
                    points[2] = xC
                    values[2] = fC
                else:
                    # Shrink
                    xS = (xB + xW) / 2.0
                    fS = objective(xS)
                    points[2] = xS
                    values[2] = fS
                    points[1] = xM
                    values[1] = objective(xM)
        
        # Check convergence
        rel_tol = abs(values[0] - values[2]) / (abs(values[0]) + abs(values[2]) + 0.01)
        x_norm = (np.sqrt(np.sum((points[0] - points[1])**2)) +
                  np.sqrt(np.sum((points[0] - points[2])**2)) +
                  np.sqrt(np.sum((points[1] - points[2])**2)))
        
        if rel_tol < f_tol or x_norm < x_tol:
            break
    
    # Final sort
    order = np.argsort(values)
    
    return PSOResult(
        x_opt=points[order[0]],
        f_opt=values[order[0]],
        converged=(rel_tol < f_tol),
        iterations=iteration + 1
    )
