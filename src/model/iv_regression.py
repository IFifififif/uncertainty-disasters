"""
Instrumental Variables Regression for BBT (2024) Model.

This module implements the two-stage least squares (2SLS) regression
used to compute moments for the GMM estimation:

First Stage:
- Regress first/second moment on disaster instruments
- Y1 = Z * beta1 (first moment regression)
- Y2 = Z * beta2 (second moment/volatility regression)

Second Stage:
- Regress GDP growth on predicted first/second moments
- growth = gamma1 * Y1_hat + gamma2 * Y2_hat + controls

Disaster instruments:
1. Natural disasters
2. Political shocks (coups)
3. Revolutions
4. Terrorist attacks

Matches MATLAB FIRST_STAGE.m and Stata ivreg2 implementation.
"""

import numpy as np
from typing import Tuple, Dict, Optional
from dataclasses import dataclass


@dataclass
class FirstStageResults:
    """Results from first stage regression."""
    
    beta_first: np.ndarray   # Coefficients for first moment
    beta_second: np.ndarray  # Coefficients for second moment
    R2_first: float
    R2_second: float
    F_stat_first: float      # F-statistic for instrument relevance
    F_stat_second: float


@dataclass  
class SecondStageResults:
    """Results from second stage regression."""
    
    gamma_first: float    # Coefficient on predicted first moment
    gamma_second: float   # Coefficient on predicted second moment
    se_first: float
    se_second: float
    R2: float


@dataclass
class IVResults:
    """Complete IV regression results."""
    
    first_stage: FirstStageResults
    second_stage: SecondStageResults
    n_obs: int


def first_stage_regression(
    Y_first: np.ndarray,
    Y_second: np.ndarray,
    Z: np.ndarray,
    X_controls: np.ndarray = None
) -> FirstStageResults:
    """
    Run first stage regression of moments on instruments.
    
    Y1 = Z * beta1 + X * delta1 + epsilon1
    Y2 = Z * beta2 + X * delta2 + epsilon2
    
    Parameters
    ----------
    Y_first : np.ndarray
        First moment (GDP growth), shape (N,).
    Y_second : np.ndarray
        Second moment (volatility), shape (N,).
    Z : np.ndarray
        Instrument matrix, shape (N, K) where K = 4 disaster types.
    X_controls : np.ndarray, optional
        Control variables, shape (N, P).
    
    Returns
    -------
    FirstStageResults
        Regression coefficients and statistics.
    """
    N = len(Y_first)
    K = Z.shape[1]
    
    # Build design matrix with controls
    if X_controls is not None:
        W = np.column_stack([np.ones(N), Z, X_controls])
    else:
        W = np.column_stack([np.ones(N), Z])
    
    # OLS for first moment
    beta_first_all = np.linalg.lstsq(W, Y_first, rcond=None)[0]
    beta_first = beta_first_all[1:K+1]  # Extract disaster coefficients
    
    # Predicted values
    Y_first_hat = W @ beta_first_all
    SS_explained = np.sum((Y_first_hat - Y_first.mean())**2)
    SS_total = np.sum((Y_first - Y_first.mean())**2)
    R2_first = SS_explained / SS_total if SS_total > 0 else 0
    
    # F-statistic for instrument relevance
    # Test H0: all disaster coefficients = 0
    Z_design = W[:, 1:K+1]
    W_restricted = W[:, 0:1]  # Just constant
    Y_first_hat_r = W_restricted @ np.linalg.lstsq(W_restricted, Y_first, rcond=None)[0]
    
    SS_r = np.sum((Y_first - Y_first_hat_r)**2)
    SS_u = np.sum((Y_first - Y_first_hat)**2)
    df_r = K
    df_denom = N - W.shape[1]
    F_stat_first = ((SS_r - SS_u) / df_r) / (SS_u / df_denom) if SS_u > 0 else 0
    
    # OLS for second moment
    beta_second_all = np.linalg.lstsq(W, Y_second, rcond=None)[0]
    beta_second = beta_second_all[1:K+1]
    
    Y_second_hat = W @ beta_second_all
    SS_explained_2 = np.sum((Y_second_hat - Y_second.mean())**2)
    SS_total_2 = np.sum((Y_second - Y_second.mean())**2)
    R2_second = SS_explained_2 / SS_total_2 if SS_total_2 > 0 else 0
    
    # F-statistic for second moment
    Y_second_hat_r = W_restricted @ np.linalg.lstsq(W_restricted, Y_second, rcond=None)[0]
    SS_r_2 = np.sum((Y_second - Y_second_hat_r)**2)
    SS_u_2 = np.sum((Y_second - Y_second_hat)**2)
    F_stat_second = ((SS_r_2 - SS_u_2) / df_r) / (SS_u_2 / df_denom) if SS_u_2 > 0 else 0
    
    return FirstStageResults(
        beta_first=beta_first,
        beta_second=beta_second,
        R2_first=R2_first,
        R2_second=R2_second,
        F_stat_first=F_stat_first,
        F_stat_second=F_stat_second
    )


def second_stage_regression(
    growth: np.ndarray,
    Y_first_hat: np.ndarray,
    Y_second_hat: np.ndarray,
    X_controls: np.ndarray = None
) -> SecondStageResults:
    """
    Run second stage regression of growth on predicted moments.
    
    growth = gamma1 * Y1_hat + gamma2 * Y2_hat + X * delta + epsilon
    
    Parameters
    ----------
    growth : np.ndarray
        GDP growth, shape (N,).
    Y_first_hat : np.ndarray
        Predicted first moment from first stage.
    Y_second_hat : np.ndarray
        Predicted second moment from first stage.
    X_controls : np.ndarray, optional
        Control variables.
    
    Returns
    -------
    SecondStageResults
        Regression coefficients and statistics.
    """
    N = len(growth)
    
    # Build design matrix
    if X_controls is not None:
        W = np.column_stack([np.ones(N), Y_first_hat, Y_second_hat, X_controls])
    else:
        W = np.column_stack([np.ones(N), Y_first_hat, Y_second_hat])
    
    # OLS
    try:
        beta = np.linalg.lstsq(W, growth, rcond=None)[0]
    except:
        beta = np.zeros(W.shape[1])
    
    gamma_first = beta[1]
    gamma_second = beta[2]
    
    # Predicted values and residuals
    growth_hat = W @ beta
    residuals = growth - growth_hat
    
    # Standard errors
    dof = N - W.shape[1]
    if dof > 0 and np.sum(residuals**2) > 0:
        sigma2 = np.sum(residuals**2) / dof
        try:
            var_beta = sigma2 * np.linalg.inv(W.T @ W)
            se_first = np.sqrt(var_beta[1, 1])
            se_second = np.sqrt(var_beta[2, 2])
        except:
            se_first = 0.0
            se_second = 0.0
    else:
        se_first = 0.0
        se_second = 0.0
    
    # R-squared
    SS_explained = np.sum((growth_hat - growth.mean())**2)
    SS_total = np.sum((growth - growth.mean())**2)
    R2 = SS_explained / SS_total if SS_total > 0 else 0
    
    return SecondStageResults(
        gamma_first=gamma_first,
        gamma_second=gamma_second,
        se_first=se_first,
        se_second=se_second,
        R2=R2
    )


def run_iv_regression(
    growth: np.ndarray,
    first_moment: np.ndarray,
    second_moment: np.ndarray,
    instruments: np.ndarray,
    controls: np.ndarray = None
) -> IVResults:
    """
    Run complete 2SLS IV regression.
    
    Parameters
    ----------
    growth : np.ndarray
        GDP growth, shape (N,).
    first_moment : np.ndarray
        First moment variable (instrumented).
    second_moment : np.ndarray
        Second moment variable (instrumented).
    instruments : np.ndarray
        Disaster instruments, shape (N, 4).
    controls : np.ndarray, optional
        Additional control variables.
    
    Returns
    -------
    IVResults
        Complete IV regression results.
    """
    N = len(growth)
    
    # First stage
    first_stage = first_stage_regression(
        first_moment, second_moment, instruments, controls
    )
    
    # Predicted moments (fitted values from first stage)
    if controls is not None:
        W = np.column_stack([np.ones(N), instruments, controls])
    else:
        W = np.column_stack([np.ones(N), instruments])
    
    # Re-compute with full coefficient vectors
    beta1_full = np.linalg.lstsq(W, first_moment, rcond=None)[0]
    beta2_full = np.linalg.lstsq(W, second_moment, rcond=None)[0]
    
    Y_first_hat = W @ beta1_full
    Y_second_hat = W @ beta2_full
    
    # Second stage
    second_stage = second_stage_regression(
        growth, Y_first_hat, Y_second_hat, controls
    )
    
    return IVResults(
        first_stage=first_stage,
        second_stage=second_stage,
        n_obs=N
    )


def compute_iv_moments(
    growth: np.ndarray,
    first_moment: np.ndarray,
    second_moment: np.ndarray,
    instruments: np.ndarray,
    sample: str = 'macro'
) -> np.ndarray:
    """
    Compute IV moments for GMM matching.
    
    Returns 10 moments:
    - 4 first-stage coefficients for first moment
    - 4 first-stage coefficients for second moment
    - 2 second-stage coefficients
    
    Parameters
    ----------
    growth : np.ndarray
        GDP growth.
    first_moment : np.ndarray
        First moment.
    second_moment : np.ndarray
        Second moment.
    instruments : np.ndarray
        Disaster instruments (4 types).
    sample : str
        'macro' or 'micro' sample.
    
    Returns
    -------
    moments : np.ndarray
        10 moments for GMM.
    """
    iv_results = run_iv_regression(
        growth, first_moment, second_moment, instruments
    )
    
    moments = np.zeros(10)
    moments[0:4] = iv_results.first_stage.beta_first
    moments[4:8] = iv_results.first_stage.beta_second
    moments[8] = iv_results.second_stage.gamma_first
    moments[9] = iv_results.second_stage.gamma_second
    
    return moments


def simulate_iv_data(
    T: int,
    disaster_probs: np.ndarray,
    disaster_levels: np.ndarray,
    disaster_unc_probs: np.ndarray,
    seed: int = 2501
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Simulate data for IV regression.
    
    This generates:
    - GDP growth series
    - First moment series (levels)
    - Second moment series (volatility)
    - Disaster indicators
    
    Parameters
    ----------
    T : int
        Number of periods.
    disaster_probs : np.ndarray
        Probability of each disaster type (4 values).
    disaster_levels : np.ndarray
        Impact on levels (4 values).
    disaster_unc_probs : np.ndarray
        Probability of uncertainty spike from disaster (4 values).
    seed : int
        Random seed.
    
    Returns
    -------
    growth : np.ndarray
    first_moment : np.ndarray
    second_moment : np.ndarray
    instruments : np.ndarray
    """
    np.random.seed(seed)
    
    # Initialize
    growth = np.zeros(T)
    first_moment = np.zeros(T)
    second_moment = np.zeros(T)
    instruments = np.zeros((T, 4))
    
    # Base parameters
    base_growth = 0.02
    base_vol = 0.02
    
    for t in range(T):
        # Check for disasters
        for d in range(4):
            if np.random.random() < disaster_probs[d]:
                instruments[t, d] = 1
                
                # Impact on growth
                growth[t] += disaster_levels[d]
                
                # Impact on volatility
                if np.random.random() < disaster_unc_probs[d]:
                    second_moment[t] += 1  # High uncertainty indicator
        
        # Base growth
        growth[t] += base_growth + np.random.normal(0, 0.01)
        
        # First moment (normalized)
        first_moment[t] = growth[t] / np.std(growth)
        
        # Second moment (log volatility)
        second_moment[t] = np.log(base_vol + second_moment[t] * 0.02)
    
    # Normalize
    first_moment = first_moment / np.std(first_moment, ddof=1)
    second_moment = second_moment / np.std(second_moment, ddof=1)
    
    return growth, first_moment, second_moment, instruments
