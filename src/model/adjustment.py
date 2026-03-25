"""
NON-CONVEX Adjustment Costs for BBT (2024) Model.

This module implements the key feature distinguishing this model from
standard DSGE models: non-convex adjustment costs for both capital and labor.

Key Features (matching Fortran VOL_GROWTH_wrapper.f90 exactly):
1. Capital partial irreversibility: Only negative investment (selling capital) is costly
2. Capital fixed adjustment cost: Lump-sum cost proportional to output when adjusting
3. Hiring costs: Linear cost when expanding labor
4. Firing costs: Linear cost when reducing labor
5. Labor fixed adjustment cost: Lump-sum cost proportional to output when adjusting labor

Original Fortran functions:
- ACk(z,a,k,l,...): Capital adjustment cost
- ACl(z,a,k,l,l_{-1},...): Labor adjustment cost
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
    deltak: float = 0.026,
    Y: float = 1.0
) -> float:
    """
    NON-CONVEX capital adjustment cost.
    
    Matches Fortran ACk function exactly.
    
    The cost includes:
    1. Partial irreversibility: Only charged when investment is NEGATIVE
       - If I < 0: ACk += |I| * capirrev
       - This captures the cost of selling capital (resale value < purchase price)
    
    2. Fixed adjustment cost: Charged when ANY non-zero investment
       - If |I| > 0: ACk += Y * capfix
       - This captures disruption costs proportional to output
    
    Fortran code (VOL_GROWTH_wrapper.f90):
    ```fortran
    ival = kprimeval - (1.0-deltak) * kval
    ACk = 0.0
    
    ! Partial irreversibility costs (only for negative investment)
    if (ival<-changetol) then
        ACk = ACk - ival * capirrev
    end if
    
    ! Fixed disruption costs (for any non-zero investment)
    if (abs(ival)>changetol) then
        yval = y(zval,aval,kval,lval,alpha,nu)
        ACk = ACk + yval * capfix
    end if
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
        Fixed adjustment cost (proportional to output).
    deltak : float
        Depreciation rate.
    Y : float
        Current output (for fixed cost calculation).
    
    Returns
    -------
    ACk : float
        Total capital adjustment cost.
    """
    changetol = 1e-10
    
    # Investment
    I = k_prime - (1 - deltak) * k
    
    # Start at 0
    ACk = 0.0
    
    # Partial irreversibility cost: only for NEGATIVE investment
    # This is the key difference from standard adjustment costs!
    if I < -changetol:
        ACk = ACk - I * capirrev  # -I = |I|
    
    # Fixed adjustment cost: for ANY non-zero investment
    if abs(I) > changetol:
        ACk = ACk + Y * capfix
    
    return ACk


def labor_adjustment_cost(
    l: float, l_prev: float, w: float,
    hirelin: float = 0.072, firelin: float = 0.072,
    labfix: float = 0.096,
    deltan: float = 0.088,
    Y: float = 1.0
) -> float:
    """
    NON-CONVEX labor adjustment cost.
    
    Matches Fortran ACl function exactly.
    
    The cost includes:
    1. Fixed adjustment cost: Only charged when labor changes
       - If |h| > 0: ACl += Y * labfix
    
    2. Hiring cost: Linear cost when expanding labor
       - If h > 0: ACl += h * w * hirelin
    
    3. Firing cost: Linear cost when reducing labor
       - If h < 0: ACl += |h| * w * firelin
    
    Fortran code (VOL_GROWTH_wrapper.f90):
    ```fortran
    hval = lval - (1.0-deltan)*lmin1val
    ACl = 0.0
    
    if (abs(hval)>changetol) then
        yval = y(zval,aval,kval,lval,alpha,nu)
        ACl = ACl + labfix * yval     ! Fixed costs
        if (hval<-changetol) ACl = ACl - hval * wval * firelin  ! Firing costs
        if (hval>changetol) ACl = ACl + hval * wval * hirelin   ! Hiring costs
    end if
    ```
    
    Parameters
    ----------
    l : float
        Current period labor.
    l_prev : float
        Previous period labor (l_{-1}).
    w : float
        Wage rate.
    hirelin : float
        Hiring cost parameter (linear).
    firelin : float
        Firing cost parameter (linear).
    labfix : float
        Fixed labor adjustment cost (proportional to output).
    deltan : float
        Labor depreciation (worker separations).
    Y : float
        Current output (for fixed cost calculation).
    
    Returns
    -------
    ACl : float
        Total labor adjustment cost.
    """
    changetol = 1e-10
    
    # Net hiring = l - (1-deltan)*l_prev
    h = l - (1 - deltan) * l_prev
    
    # Start at 0
    ACl = 0.0
    
    # Only charge costs if labor actually changes
    if abs(h) > changetol:
        # Fixed adjustment cost (proportional to output)
        ACl = ACl + labfix * Y
        
        # Firing cost (h < 0 means firing)
        if h < -changetol:
            ACl = ACl - h * w * firelin  # -h = |h|
        
        # Hiring cost (h > 0 means hiring)
        if h > changetol:
            ACl = ACl + h * w * hirelin
    
    return ACl


def compute_adjustment_costs_grid(
    k_grid: np.ndarray, l_grid: np.ndarray, w: float,
    params: ModelParameters,
    Y: float = 1.0
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
    Y : float
        Output level (for fixed costs).
    
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
                params.capirrev, params.capfix, params.deltak,
                Y=Y
            )
    
    # Labor adjustment costs
    ACl_mat = np.zeros((lnum, lnum))
    for i in range(lnum):
        for j in range(lnum):
            ACl_mat[i, j] = labor_adjustment_cost(
                l_grid[j], l_grid[i], w,
                params.hirelin, params.firelin, params.labfix, params.deltan,
                Y=Y
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
        
        # Capital adjustment cost matrix (base, without output for fixed cost)
        # Note: Fixed cost depends on Y, so we compute base cost separately
        self.ACk_base = np.zeros((p.knum, p.knum))
        for i in range(p.knum):
            for j in range(p.knum):
                # Base irreversibility cost (no output dependence)
                I = self.I_mat[i, j]
                if I < -1e-10:
                    self.ACk_base[i, j] = -I * p.capirrev  # |I| * capirrev
    
    def get_capital_adj_cost(self, k_idx: int, k_prime_idx: int, Y: float = 1.0) -> float:
        """
        Get capital adjustment cost for given indices and output.
        
        Parameters
        ----------
        k_idx : int
            Current capital index.
        k_prime_idx : int
            Next period capital index.
        Y : float
            Current output.
        
        Returns
        -------
        ACk : float
            Capital adjustment cost.
        """
        p = self.params
        
        I = self.I_mat[k_idx, k_prime_idx]
        ACk = 0.0
        
        # Partial irreversibility
        if I < -1e-10:
            ACk = ACk - I * p.capirrev
        
        # Fixed cost
        if abs(I) > 1e-10:
            ACk = ACk + Y * p.capfix
        
        return ACk
    
    def get_labor_adj_cost(
        self, l_idx: int, l_prev_idx: int, w: float, Y: float = 1.0
    ) -> float:
        """
        Get labor adjustment cost for given indices.
        
        Parameters
        ----------
        l_idx : int
            Current labor index.
        l_prev_idx : int
            Previous labor index.
        w : float
            Wage rate.
        Y : float
            Current output.
        
        Returns
        -------
        ACl : float
            Labor adjustment cost.
        """
        p = self.params
        g = self.grids
        
        return labor_adjustment_cost(
            g.l_grid[l_idx], g.l_grid[l_prev_idx], w,
            p.hirelin, p.firelin, p.labfix, p.deltan,
            Y=Y
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
        ACk = self.get_capital_adj_cost(k_idx, k_prime_idx, Y)
        
        # Labor adjustment cost
        ACl = self.get_labor_adj_cost(l_prime_idx, l_prev_idx, w, Y)
        
        # Investment
        I = self.I_mat[k_idx, k_prime_idx]
        
        # Wage bill
        WL = w * self.grids.l_grid[l_prime_idx]
        
        return price * (Y - ACk - ACl - I - WL)
