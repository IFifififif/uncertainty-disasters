"""
Shared utility functions for the BBT replication project.

Provides helpers for:
- Fixed effects demeaning (matching Stata's areg/reghdfe)
- IV/2SLS estimation (matching Stata's ivreg2)
- Output formatting (matching Stata's outreg2)
"""

import numpy as np
import pandas as pd
from typing import Optional


def demean_by_group(df: pd.DataFrame, cols: list, group_col: str) -> pd.DataFrame:
    """
    Demean columns by group (equivalent to absorbing fixed effects).

    This matches Stata's areg behavior: subtract group means from each column.
    For multiple FE, demean iteratively (Frisch-Waugh-Lovell).

    Parameters
    ----------
    df : DataFrame
    cols : list of column names to demean
    group_col : column name for grouping

    Returns
    -------
    DataFrame with demeaned columns added as '{col}_dm'
    """
    result = df.copy()
    for col in cols:
        group_means = result.groupby(group_col)[col].transform('mean')
        result[f'{col}_dm'] = result[col] - group_means
    return result


def demean_multiple_fe(
    df: pd.DataFrame,
    y_col: str,
    x_cols: list,
    fe_cols: list,
    sample_col: Optional[str] = None,
) -> tuple:
    """
    Iteratively demean by multiple fixed effects (Frisch-Waugh-Lovell).

    Matches Stata's reghdfe with absorb(fe1 fe2).

    Parameters
    ----------
    df : DataFrame
    y_col : dependent variable column
    x_cols : independent variable columns
    fe_cols : list of fixed effect column names
    sample_col : optional column to filter sample

    Returns
    -------
    y_dm : demeaned y (Series)
    X_dm : demeaned X (DataFrame)
    sample_mask : boolean mask of valid observations
    """
    data = df.copy()

    # Build sample mask
    all_cols = [y_col] + x_cols + fe_cols
    if sample_col is not None:
        data = data[data[sample_col] == 1]
    sample_mask = data[all_cols].notna().all(axis=1)
    data = data[sample_mask].copy()

    # Iterative demeaning for multiple FE
    y = data[y_col].values.copy().astype(np.float64)
    X = data[x_cols].values.copy().astype(np.float64)

    for fe_col in fe_cols:
        groups = data[fe_col].values
        unique_groups = np.unique(groups)

        for g in unique_groups:
            mask = groups == g
            n_g = mask.sum()
            if n_g > 0:
                y_mean = y[mask].mean()
                x_mean = X[mask].mean(axis=0)
                y[mask] -= y_mean
                X[mask] -= x_mean

    return pd.Series(y, index=data.index), pd.DataFrame(X, columns=x_cols, index=data.index), sample_mask


def ols_with_cluster_se(
    y: np.ndarray,
    X: np.ndarray,
    clusters: np.ndarray,
) -> dict:
    """
    OLS estimation with cluster-robust standard errors.

    Matches Stata's areg ..., cluster(var).

    Parameters
    ----------
    y : (N,) dependent variable
    X : (N, K) regressors (should NOT include constant if demeaned)
    clusters : (N,) cluster identifiers

    Returns
    -------
    dict with keys: coef, se, tstat, pval, nobs, nclusters, r2, residuals
    """
    # OLS: beta = (X'X)^{-1} X'y
    XtX = X.T @ X
    Xty = X.T @ y
    try:
        XtX_inv = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        XtX_inv = np.linalg.pinv(XtX)

    beta = XtX_inv @ Xty
    residuals = y - X @ beta
    N, K = X.shape

    # Cluster-robust SEs
    unique_clusters = np.unique(clusters)
    n_clusters = len(unique_clusters)

    # Meat: sum of (X_g' e_g)(X_g' e_g)' / (1 - K/N)^2 * N/(N-1) * nc/(nc-1)
    meat = np.zeros((K, K))
    for c in unique_clusters:
        mask = clusters == c
        Xc = X[mask]
        ec = residuals[mask]
        Xu_e = Xc.T @ ec
        meat += np.outer(Xu_e, Xu_e)

    # Small sample correction
    bread = XtX_inv
    V = bread @ meat @ bread

    # Degrees of freedom correction
    df_correction = (N / (N - 1)) * (n_clusters / (n_clusters - 1))
    V *= df_correction

    se = np.sqrt(np.diag(V))

    # R-squared
    y_mean = y.mean()
    ss_tot = np.sum((y - y_mean) ** 2)
    ss_res = np.sum(residuals ** 2)
    r2 = 1 - ss_res / ss_tot

    return {
        'coef': beta,
        'se': se,
        'residuals': residuals,
        'nobs': N,
        'nclusters': n_clusters,
        'r2': r2,
        'V': V,
    }


def iv2sls_with_cluster_se(
    y: np.ndarray,
    X_endog: np.ndarray,
    X_exog: np.ndarray,
    Z: np.ndarray,
    clusters: np.ndarray,
    partial_out: Optional[np.ndarray] = None,
) -> dict:
    """
    Two-Stage Least Squares with cluster-robust standard errors.

    Matches Stata's ivreg2 with cluster() and partial() options.

    Parameters
    ----------
    y : (N,) dependent variable
    X_endog : (N, L_endog) endogenous regressors
    X_exog : (N, L_exog) exogenous regressors (included in both stages)
    Z : (N, L_z) excluded instruments
    clusters : (N,) cluster identifiers
    partial_out : (N, L_partial) variables to partial out (FWL)

    Returns
    -------
    dict with coef, se, first_stage results, Hansen J, etc.
    """
    N = len(y)

    # Step 1: Partial out if specified (FWL theorem)
    if partial_out is not None:
        # Use iterative demeaning for each column of partial_out
        # This handles collinear dummy variables (e.g., cc* + yy*)
        # by demeaning one group at a time (Frisch-Waugh-Lovell)
        y_res = y.copy()
        X_endog_res = X_endog.copy()
        X_exog_res = X_exog.copy()
        Z_res = Z.copy()

        # partial_out columns are grouped: first cc_cols, then yy_cols
        # We need to identify group boundaries from the data
        # Since partial_out is [cc* yy*], we demean by each group
        # But we don't have group labels here — use the column structure
        # Each column is a dummy; we identify groups by finding columns
        # that are mutually exclusive (exactly one is 1 per group)

        # Alternative: use QR-based projection which handles rank deficiency
        P = partial_out
        Q, R = np.linalg.qr(P, mode='reduced')
        # Determine numerical rank
        diag_R = np.abs(np.diag(R))
        rank = np.sum(diag_R > 1e-10 * diag_R[0])
        Q_rank = Q[:, :rank]

        # Project out using full-rank QR decomposition
        Px = Q_rank @ Q_rank.T
        y_res = y - Px @ y
        if X_endog_res.ndim == 1:
            X_endog_res = X_endog_res - Px @ X_endog_res
        else:
            X_endog_res = X_endog_res - Px @ X_endog_res
        if X_exog_res.shape[1] > 0:
            X_exog_res = X_exog_res - Px @ X_exog_res
        Z_res = Z - Px @ Z
    else:
        y_res = y.copy()
        X_endog_res = X_endog.copy()
        X_exog_res = X_exog.copy()
        Z_res = Z.copy()

    # Build full instrument matrix: [X_exog, Z]
    W = np.hstack([X_exog_res, Z_res])

    # First stage: regress each endogenous variable on W
    L_endog = X_endog_res.shape[1] if X_endog_res.ndim > 1 else 1
    if X_endog_res.ndim == 1:
        X_endog_res = X_endog_res.reshape(-1, 1)

    first_stage_results = []
    X_hat_list = []

    # Check for and remove constant/near-constant columns in W
    W_stds = np.std(W, axis=0)
    valid_cols = W_stds > 1e-10
    if not np.all(valid_cols):
        # Keep only non-constant columns
        W = W[:, valid_cols]
        removed_count = np.sum(~valid_cols)
        if removed_count > 0:
            import warnings
            warnings.warn(f"Removed {removed_count} constant instrument(s) after partialling out FE")

    for j in range(L_endog):
        # OLS of X_endog_j on W using pseudo-inverse for numerical stability
        WtW = W.T @ W
        try:
            WtW_inv = np.linalg.inv(WtW)
        except np.linalg.LinAlgError:
            # Use pseudo-inverse if matrix is singular
            WtW_inv = np.linalg.pinv(WtW)
        pi_j = WtW_inv @ (W.T @ X_endog_res[:, j])
        X_hat_j = W @ pi_j
        X_hat_list.append(X_hat_j)

        # First stage F-stat
        resid_j = X_endog_res[:, j] - X_hat_j
        ssr_restricted = np.sum(resid_j ** 2)
        # Unrestricted: regress on X_exog only
        if X_exog_res.shape[1] > 0:
            Xe_Xe_inv = np.linalg.inv(X_exog_res.T @ X_exog_res)
            pi_r = Xe_Xe_inv @ (X_exog_res.T @ X_endog_res[:, j])
            resid_r = X_endog_res[:, j] - X_exog_res @ pi_r
            ssr_unrestricted = np.sum(resid_r ** 2)
        else:
            ssr_unrestricted = np.sum(X_endog_res[:, j] ** 2)

        n_instr = Z_res.shape[1]
        F_stat = ((ssr_unrestricted - ssr_restricted) / n_instr) / (ssr_restricted / (N - W.shape[1]))

        first_stage_results.append({
            'coef': pi_j,
            'F_stat': F_stat,
            'X_hat': X_hat_j,
        })

    X_hat = np.column_stack(X_hat_list)

    # Second stage: regress y on [X_exog, X_hat]
    if X_exog_res.shape[1] == 0:
        X_2sls = X_hat
    else:
        X_2sls = np.hstack([X_exog_res, X_hat])
    result = ols_with_cluster_se(y_res, X_2sls, clusters)

    # Organize coefficients
    n_exog = X_exog_res.shape[1]
    coef_endog = result['coef'][n_exog:]
    se_endog = result['se'][n_exog:]

    # Hansen J statistic (overidentification test)
    # J = N * e' P_Z e / (e'e/N)
    e = result['residuals']
    try:
        Pz = W @ np.linalg.inv(W.T @ W) @ W.T
    except np.linalg.LinAlgError:
        Pz = W @ np.linalg.pinv(W.T @ W) @ W.T
    e_Pz_e = e.T @ Pz @ e
    e_e = e.T @ e
    J_stat = N * e_Pz_e / (e_e / N)
    # Use the actual number of instruments used (after removing constants)
    n_instr_used = W.shape[1]
    n_overid = n_instr_used - L_endog
    if n_overid > 0:
        from scipy.stats import chi2
        J_pval = 1 - chi2.cdf(J_stat, n_overid)
    else:
        J_pval = np.nan

    return {
        'coef': result['coef'],
        'se': result['se'],
        'coef_endog': coef_endog,
        'se_endog': se_endog,
        'first_stage': first_stage_results,
        'J_stat': J_stat,
        'J_pval': J_pval,
        'nobs': N,
        'nclusters': len(np.unique(clusters)),
        'residuals': result['residuals'],
    }


def get_cc_yy_cols(df: pd.DataFrame) -> tuple:
    """
    Extract country dummies (cc*) and time dummies (yy*) from the dataframe.

    Returns
    -------
    cc_cols : list of country dummy column names
    yy_cols : list of year-quarter dummy column names
    """
    cc_cols = sorted([c for c in df.columns if c.startswith('cc')])
    yy_cols = sorted([c for c in df.columns if c.startswith('yy')],
                     key=lambda x: int(x[2:]))
    return cc_cols, yy_cols


def format_coef_table(
    coefs: np.ndarray,
    ses: np.ndarray,
    var_names: list,
    title: str = "",
    nobs: int = 0,
    nclusters: int = 0,
    addtext: Optional[dict] = None,
) -> str:
    """
    Format regression results as a string table (matching outreg2 style).

    Parameters
    ----------
    coefs : coefficient array
    ses : standard error array
    var_names : variable names
    title : table title
    nobs : number of observations
    nclusters : number of clusters
    addtext : dict of additional text entries

    Returns
    -------
    Formatted string
    """
    lines = []
    if title:
        lines.append(title)
        lines.append("=" * 60)

    lines.append(f"{'Variable':<30} {'Coef':>10} {'Std.Err':>10} {'t-stat':>10}")
    lines.append("-" * 60)

    for i, name in enumerate(var_names):
        if i < len(coefs):
            tstat = coefs[i] / ses[i] if ses[i] > 0 else 0
            sig = ""
            if abs(tstat) > 2.576:
                sig = "***"
            elif abs(tstat) > 1.960:
                sig = "**"
            elif abs(tstat) > 1.645:
                sig = "*"
            lines.append(f"{name:<30} {coefs[i]:>10.4f} {ses[i]:>10.4f} {tstat:>10.4f}{sig}")

    lines.append("-" * 60)
    if nobs:
        lines.append(f"{'Observations':<30} {nobs:>10d}")
    if nclusters:
        lines.append(f"{'Clusters':<30} {nclusters:>10d}")

    if addtext:
        for k, v in addtext.items():
            lines.append(f"{k:<30} {v:>10}")

    return "\n".join(lines)
