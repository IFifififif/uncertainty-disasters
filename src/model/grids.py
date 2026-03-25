"""
State Space Grid Construction for BBT (2024) Model.

Implements the 5-dimensional state space:
- z: Idiosyncratic productivity (Tauchen discretization)
- a: Aggregate productivity (Tauchen discretization)
- s: Volatility state (2 states: low/high uncertainty)
- k: Capital stock
- l: Labor stock

Plus forecast states for aggregate capital (Kbar).

Matches Fortran VOL_GROWTH_wrapper.f90 exactly.
"""

import numpy as np
from scipy.stats import norm
from typing import Tuple, Dict
from dataclasses import dataclass

from .params import ModelParameters


@dataclass
class StateGrids:
    """Container for all state space grids and transition matrices."""
    
    # Grids
    k_grid: np.ndarray          # Capital grid
    l_grid: np.ndarray          # Labor grid
    z_grid: np.ndarray          # Idiosyncratic productivity grid
    a_grid: np.ndarray          # Aggregate productivity grid
    kbar_grid: np.ndarray       # Aggregate capital forecast grid
    
    # Volatility levels
    sigmaz_grid: np.ndarray     # Idio vol: [low, high]
    sigmaa_grid: np.ndarray     # Agg vol: [low, high]
    
    # Transition matrices
    pr_mat_z: np.ndarray        # Shape: (znum, znum, snum)
    pr_mat_a: np.ndarray        # Shape: (anum, anum, snum)
    pr_mat_s: np.ndarray        # Shape: (snum, snum)
    pr_mat_full: np.ndarray     # Combined transition matrix
    
    # Indexing arrays (matching Fortran)
    exog_pos: np.ndarray        # Exogenous state positions, shape (numexog, 4)
    endog_pos: np.ndarray       # Endogenous state positions, shape (numendog, 2)
    
    # Dimensions
    numexog: int
    numendog: int
    numstates: int
    numexog_next: int           # Number of next period exogenous states


def build_grids(params: ModelParameters) -> StateGrids:
    """
    Build all state space grids and transition matrices.
    
    Matches Fortran VOL_GROWTH_wrapper.f90 exactly.
    
    Parameters
    ----------
    params : ModelParameters
        Model parameters instance.
    
    Returns
    -------
    StateGrids
        Container with all grids and transition matrices.
    """
    p = params
    
    # =====================
    # Build Capital Grid
    # =====================
    # Exponential spacing respecting depreciation
    # Fortran line 380: kmax = exp(log(kmin) - dble(knum-1) * log(1.0-deltak))
    kmax = np.exp(np.log(p.kmin) - (p.knum - 1) * np.log(1 - p.deltak))
    # Fortran line 500: call linspace(k0,log(kmin),log(kmax),knum); k0=exp(k0)
    k_grid = np.exp(np.linspace(np.log(p.kmin), np.log(kmax), p.knum))
    
    # =====================
    # Build Labor Grid
    # =====================
    # Fortran line 384: lmax = exp(log(lmin) - dble(lnum-1) * log(1.0-deltan))
    lmax = np.exp(np.log(p.lmin) - (p.lnum - 1) * np.log(1 - p.deltan))
    l_grid = np.exp(np.linspace(np.log(p.lmin), np.log(lmax), p.lnum))
    
    # =====================
    # Build Aggregate Capital Forecast Grid
    # =====================
    # Fortran line 504: kbar0 = (/kbarmin,kbarmax/)
    kbar_grid = np.array([p.kbarmin, p.kbarmax])
    
    # =====================
    # Build Volatility Grids
    # =====================
    # Fortran lines 507-508
    sigmaz_grid = np.array([p.sigmaz, p.zjump * p.sigmaz])
    sigmaa_grid = np.array([p.sigmaa, p.ajump * p.sigmaa])
    
    # =====================
    # Build Productivity Grids (Tauchen method)
    # =====================
    # Idiosyncratic productivity z: log-space AR(1)
    # Fortran line 519: call unctauchen(pr_mat_z,z0,znum,zmin,zmax,snum,sigmazgrid,rhoz,dble(0.0))
    z_grid, pr_mat_z = build_tauchen_grid_with_sv(
        p.zmin, p.zmax, p.znum, p.rhoz, sigmaz_grid, p.snum,
        meanshift=0.0
    )
    
    # Aggregate productivity a: level-space AR(1)
    # Fortran line 520: call unctauchen(pr_mat_a,a0,anum,amin,amax,snum,sigmaagrid,rhoa,meanshifta)
    a_grid = np.linspace(p.amin, p.amax, p.anum)
    pr_mat_a = build_tauchen_transition_sv(
        a_grid, p.anum, p.rhoa, sigmaa_grid, p.snum,
        meanshift=0.0  # Base case; adjusted during simulation
    )
    
    # =====================
    # Build Uncertainty Transition Matrix
    # =====================
    # Fortran lines 509-510
    pr_mat_s = np.zeros((p.snum, p.snum))
    pr_mat_s[0, :] = [1 - p.uncfreq, p.uncfreq]
    pr_mat_s[0, :] /= pr_mat_s[0, :].sum()
    pr_mat_s[1, :] = [1 - p.uncpers, p.uncpers]
    pr_mat_s[1, :] /= pr_mat_s[1, :].sum()
    
    # =====================
    # Build Full State Indexing
    # =====================
    # Exogenous states: (z, a, s, s_{-1})
    # Fortran line 388: numexog = znum*anum*snum*snum
    numexog = p.znum * p.anum * p.snum * p.snum
    exog_pos = np.zeros((numexog, 4), dtype=int)
    
    # Fortran lines 523-543
    ct = 0
    for zct in range(p.znum):
        for act in range(p.anum):
            for sct in range(p.snum):
                for sm1ct in range(p.snum):
                    exog_pos[ct, 0] = zct
                    exog_pos[ct, 1] = act
                    exog_pos[ct, 2] = sct
                    exog_pos[ct, 3] = sm1ct
                    ct += 1
    
    # Endogenous states: (k, l_{-1})
    # Fortran line 389: numendog = knum*lnum
    numendog = p.knum * p.lnum
    endog_pos = np.zeros((numendog, 2), dtype=int)
    
    # Fortran lines 581-586
    ct = 0
    for kct in range(p.knum):
        for lmin1ct in range(p.lnum):
            endog_pos[ct, 0] = kct
            endog_pos[ct, 1] = lmin1ct
            ct += 1
    
    # =====================
    # Build Full Transition Matrix
    # =====================
    # Fortran lines 534-541
    pr_mat_full, numexog_next = build_full_transition_matrix(
        pr_mat_z, pr_mat_a, pr_mat_s, p.znum, p.anum, p.snum
    )
    
    # Fortran line 391: numstates = numexog * numendog * numfcst
    numstates = numexog * numendog * p.kbarnum
    
    return StateGrids(
        k_grid=k_grid,
        l_grid=l_grid,
        z_grid=z_grid,
        a_grid=a_grid,
        kbar_grid=kbar_grid,
        sigmaz_grid=sigmaz_grid,
        sigmaa_grid=sigmaa_grid,
        pr_mat_z=pr_mat_z,
        pr_mat_a=pr_mat_a,
        pr_mat_s=pr_mat_s,
        pr_mat_full=pr_mat_full,
        exog_pos=exog_pos,
        endog_pos=endog_pos,
        numexog=numexog,
        numendog=numendog,
        numstates=numstates,
        numexog_next=numexog_next
    )


def build_tauchen_grid_with_sv(
    zmin: float, zmax: float, znum: int, 
    rho: float, sigma_grid: np.ndarray, snum: int,
    meanshift: float = 0.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build productivity grid using Tauchen method with stochastic volatility.
    
    Matches Fortran `unctauchen` subroutine exactly.
    
    Parameters
    ----------
    zmin, zmax : float
        Grid bounds.
    znum : int
        Number of grid points.
    rho : float
        Autocorrelation.
    sigma_grid : np.ndarray
        Volatility levels for each state (unconditional std in log space).
    snum : int
        Number of volatility states.
    meanshift : float
        Mean shift adjustment (for aggregate productivity only).
    
    Returns
    -------
    z_grid : np.ndarray
        Productivity grid.
    pr_mat : np.ndarray
        Transition matrix, shape (znum, znum, snum).
    """
    # Build grid in log space (Fortran: call linspace(log_grid,log(zmin),log(zmax),znum))
    log_grid = np.linspace(np.log(zmin), np.log(zmax), znum)
    z_grid = np.exp(log_grid)
    
    # Grid spacing (Fortran: gridinc = log_grid(2)-log_grid(1))
    gridinc = log_grid[1] - log_grid[0]
    
    # Build transition matrix for each volatility state
    pr_mat = np.zeros((znum, znum, snum))
    
    for s in range(snum):
        # Fortran uses sigma_grid directly as the std dev in normcdf
        standdev = sigma_grid[s]
        
        for i in range(znum):
            # Conditional mean: rho * log_grid[i] + meanshift
            mu = rho * log_grid[i] + meanshift
            
            # Middle intervals (zprimect = 2 to znum-1)
            for j in range(1, znum - 1):
                pr_mat[i, j, s] = (
                    norm.cdf(log_grid[j] + gridinc / 2, mu, standdev) -
                    norm.cdf(log_grid[j] - gridinc / 2, mu, standdev)
                )
            
            # First interval (Fortran: transarray(zct,1,sct))
            pr_mat[i, 0, s] = norm.cdf(log_grid[0] + gridinc / 2, mu, standdev)
            
            # Last interval (Fortran: transarray(zct,znum,sct))
            pr_mat[i, znum - 1, s] = 1 - norm.cdf(log_grid[znum - 1] - gridinc / 2, mu, standdev)
            
            # Normalize to sum to 1
            pr_mat[i, :, s] /= pr_mat[i, :, s].sum()
    
    return z_grid, pr_mat


def build_tauchen_transition_sv(
    a_grid: np.ndarray, anum: int, rho: float,
    sigma_grid: np.ndarray, snum: int,
    meanshift: float = 0.0
) -> np.ndarray:
    """
    Build transition matrix for aggregate productivity with stochastic volatility.
    
    Matches Fortran `unctauchen` subroutine for level-space grid.
    
    Parameters
    ----------
    a_grid : np.ndarray
        Productivity grid (in level space, not log space).
    anum : int
        Number of grid points.
    rho : float
        Autocorrelation.
    sigma_grid : np.ndarray
        Volatility levels for each state.
    snum : int
        Number of volatility states.
    meanshift : float
        Mean shift adjustment for disaster-induced mean impact.
    
    Returns
    -------
    pr_mat : np.ndarray
        Transition matrix, shape (anum, anum, snum).
    """
    pr_mat = np.zeros((anum, anum, snum))
    
    # Grid spacing (uniform in level space)
    gridinc = a_grid[1] - a_grid[0]
    
    for s in range(snum):
        # Fortran uses sigma_grid directly as the std dev
        standdev = sigma_grid[s]
        
        for i in range(anum):
            # Conditional mean: rho * a_grid[i] + meanshift
            mu = rho * a_grid[i] + meanshift
            
            # Middle intervals
            for j in range(1, anum - 1):
                pr_mat[i, j, s] = (
                    norm.cdf(a_grid[j] + gridinc / 2, mu, standdev) -
                    norm.cdf(a_grid[j] - gridinc / 2, mu, standdev)
                )
            
            # First interval
            pr_mat[i, 0, s] = norm.cdf(a_grid[0] + gridinc / 2, mu, standdev)
            
            # Last interval
            pr_mat[i, anum - 1, s] = 1 - norm.cdf(a_grid[anum - 1] - gridinc / 2, mu, standdev)
            
            # Normalize
            pr_mat[i, :, s] /= pr_mat[i, :, s].sum()
    
    return pr_mat


def build_full_transition_matrix(
    pr_mat_z: np.ndarray, pr_mat_a: np.ndarray, pr_mat_s: np.ndarray,
    znum: int, anum: int, snum: int
) -> Tuple[np.ndarray, int]:
    """
    Build combined transition matrix over (z, a, s) states.
    
    Matches Fortran lines 534-541 exactly.
    
    Note: s_{-1} (lagged uncertainty) is known when transitioning,
    so it affects indexing but not the transition probabilities.
    
    Current state: (z, a, s, s_{-1})
    Next state: (z', a', s'), where s_{-1}' = s (current period's s)
    
    Parameters
    ----------
    pr_mat_z : np.ndarray
        Idiosyncratic productivity transition, shape (znum, znum, snum).
    pr_mat_a : np.ndarray
        Aggregate productivity transition, shape (anum, anum, snum).
    pr_mat_s : np.ndarray
        Uncertainty transition, shape (snum, snum).
    znum, anum, snum : int
        Grid sizes.
    
    Returns
    -------
    pr_mat_full : np.ndarray
        Combined transition matrix, shape (numexog, numexog_next).
    numexog_next : int
        Number of next period exogenous states.
    """
    numexog = znum * anum * snum * snum
    numexog_next = znum * anum * snum  # s_{-1}' is determined by current s
    
    pr_mat_full = np.zeros((numexog, numexog_next))
    
    ct = 0
    for zct in range(znum):
        for act in range(anum):
            for sct in range(snum):
                for sm1ct in range(snum):
                    # Current state index
                    current_idx = ct
                    ct += 1
                    
                    # Next period: s becomes s_{-1}'
                    # Fortran line 536: primect=(zprimect-1)*anum*snum*snum + (aprimect-1)*snum*snum + (sprimect-1)*snum + sct
                    for zpct in range(znum):
                        for apct in range(anum):
                            for spct in range(snum):
                                # Next period index
                                next_idx = zpct * anum * snum + apct * snum + spct
                                
                                # Probability: P(z'|z,s) * P(a'|a,s) * P(s'|s)
                                # Fortran line 537
                                pr_mat_full[current_idx, next_idx] = (
                                    pr_mat_z[zct, zpct, sct] *
                                    pr_mat_a[act, apct, sct] *
                                    pr_mat_s[sct, spct]
                                )
    
    # Normalize rows (Fortran line 541)
    for i in range(numexog):
        if pr_mat_full[i, :].sum() > 0:
            pr_mat_full[i, :] /= pr_mat_full[i, :].sum()
    
    return pr_mat_full, numexog_next


def get_state_indices(
    exogct: int, endogct: int, fcstct: int,
    grids: StateGrids, params: ModelParameters
) -> Dict[str, int]:
    """
    Extract individual state indices from combined state index.
    
    Parameters
    ----------
    exogct, endogct, fcstct : int
        State indices.
    grids : StateGrids
        Grid container.
    params : ModelParameters
        Parameters.
    
    Returns
    -------
    indices : dict
        Dictionary with zct, act, sct, sm1ct, kct, lct.
    """
    return {
        'zct': grids.exog_pos[exogct, 0],
        'act': grids.exog_pos[exogct, 1],
        'sct': grids.exog_pos[exogct, 2],
        'sm1ct': grids.exog_pos[exogct, 3],
        'kct': grids.endog_pos[endogct, 0],
        'lct': grids.endog_pos[endogct, 1],
        'fcstct': fcstct
    }
