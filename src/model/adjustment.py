"""
NON-CONVEX Adjustment Costs for BBT (2024) Model.

This module implements the key feature distinguishing this model from
standard DSGE models: non-convex adjustment costs for both capital and labor.

Key Features:
1. Capital irreversibility: Selling capital is costly
2. Fixed adjustment costs: Lump-sum cost when adjusting
3. Hiring costs: Linear cost when expanding labor
4. Firing costs: Linear cost when reducing labor

Original Fortran code:
- ACk: Capital adjustment cost function
- ACl: Labor adjustment cost function
"""

import numpy as np
from typing import Union
from .params import ModelParameters


def output(z: float, a: float, k: float, l: float, 
           alpha: float = 0.25, nu: float = 0.5) -> float:
    """
    Production function: Y = z * a * k^alpha * l^nu
    
    Parameters
    ----------
    z : float
        Idiosyncratic productivity.
    a : float
        Aggregate productivity.
    k : float
        Capital.
    l : float
        Labor.
    alpha, nu : float
        Factor shares.
    
    Returns
    -------
    Y : float
        Output.
    """
    return z * a * (k ** alpha) * (l ** nu)


def capital_adjustment_cost(
    k_prime: float, k: float,
    capirrev: float = 0.339, capfix: float = 0.0,
    deltak: float = 0.026
) -> float:
    """
    NON-CONVEX capital adjustment cost.
    
    This is the CRITICAL difference from standard quadratic adjustment costs!
    
    The cost includes:
    1. Purchase cost for new capital: max(I, 0) = max(k' - (1-delta)*k, 0)
    2. Irreversibility cost: capirrev * |k' - k| when adjusting
    3. Fixed cost: capfix * k when adjusting (if capfix > 0)
    
    Fortran code (base_lib.f90):
    ```fortran
    ACk = capirrev * abs(k_new - k_old) + capfix * I(k_new != k_old)
    ```
    
    Parameters
    ----------
    k_prime : float
        Next period capital.
    k : float
        Current capital.
    capirrev : float
        Irreversibility cost parameter (partial irreversibility).
    capfix : float
        Fixed adjustment cost.
    deltak : float
        Depreciation rate.
    
    Returns
    -------
    ACk : float
        Total capital adjustment cost.
    """
    # Investment
    I = k_prime - (1 - deltak) * k
    
    # Adjustment magnitude
    adj = np.abs(k_prime - k)
    
    # Is there an adjustment? (tolerance for numerical precision)
    is_adjusting = adj > 1e-10
    
    # Irreversibility cost: linear cost of adjustment
    # This captures partial irreversibility - selling capital is costly
    ac_irrev = capirrev * adj
    
    # Fixed cost: only incurred if adjusting
    ac_fix = capfix * k * is_adjusting
    
    return ac_irrev + ac_fix


def labor_adjustment_cost(
    l_prime: float, l: float, w: float,
    hirelin: float = 0.072, firelin: float = 0.072,
    labfix: float = 0.096
) -> float:
    """
    NON-CONVEX labor adjustment cost.
    
    The cost includes:
    1. Hiring cost: hirelin * max(l' - l, 0) * w
    2. Firing cost: firelin * max(l - l', 0) * w
    3. Fixed cost: labfix * w when adjusting
    
    Fortran code (base_lib.f90):
    ```fortran
    ACl = hirelin * max(l_new - l_old, 0) + firelin * max(l_old - l_new, 0) 
        + labfix * I(l_new != l_old)
    ACl = ACl * w
    ```
    
    Parameters
    ----------
    l_prime : float
        Next period labor.
    l : float
        Current labor.
    w : float
        Wage rate.
    hirelin : float
        Hiring cost parameter (linear).
    firelin : float
        Firing cost parameter (linear).
    labfix : float
        Fixed labor adjustment cost.
    
    Returns
    -------
    ACl : float
        Total labor adjustment cost.
    """
    # Labor change
    dl = l_prime - l
    
    # Hiring cost (only if expanding)
    ac_hire = hirelin * max(dl, 0)
    
    # Firing cost (only if contracting)
    ac_fire = firelin * max(-dl, 0)
    
    # Fixed cost: only incurred if adjusting
    is_adjusting = np.abs(dl) > 1e-10
    ac_fix = labfix * is_adjusting
    
    return (ac_hire + ac_fire + ac_fix) * w


def compute_adjustment_costs_grid(
    k_grid: np.ndarray, l_grid: np.ndarray, w: float,
    params: ModelParameters
) -> tuple:
    """
    Pre-compute adjustment cost matrices for efficiency.
    
    Parameters
    ----------
    k_grid : np.ndarray
        Capital grid, shape (knum,).
    l_grid : np.ndarray
        Labor grid, shape (lnum,).
    w : float
        Wage rate.
    params : ModelParameters
        Model parameters.
    
    Returns
    -------
    ACk_mat : np.ndarray
        Capital adjustment costs, shape (knum, knum).
        ACk_mat[i, j] = cost of moving from k_grid[i] to k_grid[j].
    ACl_mat : np.ndarray
        Labor adjustment costs, shape (lnum, lnum).
        ACl_mat[i, j] = cost of moving from l_grid[i] to l_grid[j].
    """
    knum = len(k_grid)
    lnum = len(l_grid)
    
    # Capital adjustment costs
    ACk_mat = np.zeros((knum, knum))
    for i in range(knum):
        for j in range(knum):
            ACk_mat[i, j] = capital_adjustment_cost(
                k_grid[j], k_grid[i],
                params.capirrev, params.capfix, params.deltak
            )
    
    # Labor adjustment costs
    ACl_mat = np.zeros((lnum, lnum))
    for i in range(lnum):
        for j in range(lnum):
            ACl_mat[i, j] = labor_adjustment_cost(
                l_grid[j], l_grid[i], w,
                params.hirelin, params.firelin, params.labfix
            )
    
    return ACk_mat, ACl_mat


def compute_investment_matrix(
    k_grid: np.ndarray, deltak: float
) -> np.ndarray:
    """
    Pre-compute investment matrix.
    
    Parameters
    ----------
    k_grid : np.ndarray
        Capital grid.
    deltak : float
        Depreciation rate.
    
    Returns
    -------
    I_mat : np.ndarray
        Investment matrix, shape (knum, knum).
        I_mat[i, j] = k_grid[j] - (1-deltak) * k_grid[i].
    """
    knum = len(k_grid)
    I_mat = np.zeros((knum, knum))
    
    for i in range(knum):
        for j in range(knum):
            I_mat[i, j] = k_grid[j] - (1 - deltak) * k_grid[i]
    
    return I_mat


def compute_output_matrix(
    z_grid: np.ndarray, a_grid: np.ndarray,
    k_grid: np.ndarray, l_grid: np.ndarray,
    alpha: float, nu: float
) -> np.ndarray:
    """
    Pre-compute output matrix.
    
    Parameters
    ----------
    z_grid, a_grid, k_grid, l_grid : np.ndarray
        State grids.
    alpha, nu : float
        Factor shares.
    
    Returns
    -------
    Y_mat : np.ndarray
        Output matrix, shape (znum, anum, knum, lnum).
    """
    znum = len(z_grid)
    anum = len(a_grid)
    knum = len(k_grid)
    lnum = len(l_grid)
    
    Y_mat = np.zeros((znum, anum, knum, lnum))
    
    for zi in range(znum):
        for ai in range(anum):
            for ki in range(knum):
                for li in range(lnum):
                    Y_mat[zi, ai, ki, li] = output(
                        z_grid[zi], a_grid[ai],
                        k_grid[ki], l_grid[li],
                        alpha, nu
                    )
    
    return Y_mat


class AdjustmentCostCalculator:
    """
    Calculator for adjustment costs with pre-computed matrices.
    
    This class pre-computes all adjustment cost matrices for efficiency,
    matching the Fortran approach.
    """
    
    def __init__(self, params: ModelParameters, grids):
        """
        Initialize the calculator.
        
        Parameters
        ----------
        params : ModelParameters
            Model parameters.
        grids : StateGrids
            State grids container.
        """
        self.params = params
        self.grids = grids
        
        # Pre-compute matrices
        self._build_matrices()
    
    def _build_matrices(self):
        """Build all pre-computed matrices."""
        p = self.params
        g = self.grids
        
        # Investment matrix
        self.I_mat = compute_investment_matrix(g.k_grid, p.deltak)
        
        # Output matrix
        self.Y_mat = compute_output_matrix(
            g.z_grid, g.a_grid, g.k_grid, g.l_grid,
            p.alpha, p.nu
        )
        
        # Capital adjustment cost matrix (knum x knum)
        self.ACk_base = np.zeros((p.knum, p.knum))
        for i in range(p.knum):
            for j in range(p.knum):
                self.ACk_base[i, j] = capital_adjustment_cost(
                    g.k_grid[j], g.k_grid[i],
                    p.capirrev, p.capfix, p.deltak
                )
    
    def get_labor_adj_cost(self, l_prime_idx: int, l_idx: int, w: float) -> float:
        """
        Get labor adjustment cost for given indices and wage.
        
        Parameters
        ----------
        l_prime_idx : int
            Next period labor index.
        l_idx : int
            Current labor index.
        w : float
            Wage rate.
        
        Returns
        -------
        ACl : float
            Labor adjustment cost.
        """
        p = self.params
        g = self.grids
        
        return labor_adjustment_cost(
            g.l_grid[l_prime_idx], g.l_grid[l_idx], w,
            p.hirelin, p.firelin, p.labfix
        )
    
    def get_period_return(
        self, z_idx: int, a_idx: int, k_idx: int, l_idx: int,
        k_prime_idx: int, l_prime_idx: int, l_prev_idx: int,
        price: float, w: float
    ) -> float:
        """
        Compute period return (profit) for given state and policy.
        
        Return = p * (Y - ACk - ACl - I - w*l)
        
        Parameters
        ----------
        z_idx, a_idx, k_idx, l_idx : int
            State indices.
        k_prime_idx, l_prime_idx : int
            Policy indices.
        l_prev_idx : int
            Lagged labor index (for labor adjustment cost).
        price : float
            Output price.
        w : float
            Wage rate.
        
        Returns
        -------
        ret : float
            Period return.
        """
        p = self.params
        
        # Output
        Y = self.Y_mat[z_idx, a_idx, k_idx, l_idx]
        
        # Capital adjustment cost
        ACk = self.ACk_base[k_idx, k_prime_idx]
        
        # Labor adjustment cost
        ACl = self.get_labor_adj_cost(l_prime_idx, l_prev_idx, w)
        
        # Investment
        I = self.I_mat[k_idx, k_prime_idx]
        
        # Wage bill
        WL = w * self.grids.l_grid[l_prime_idx]
        
        return price * (Y - ACk - ACl - I - WL)
