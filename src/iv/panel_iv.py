"""
Panel IV Regressions for BBT (2024).

Replicates Stata code: Panel IV Code.do
Produces Tables 1-6 from the paper.

Original: Stata/MP 15.1, ~1 min runtime
Python: Uses linearmodels + manual FE demeaning for bit-for-bit match
"""

import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.utils.regression import (
    demean_multiple_fe,
    iv2sls_with_cluster_se,
    ols_with_cluster_se,
    get_cc_yy_cols,
    format_coef_table,
)


class PanelIV:
    """
    Replicates the Panel IV regressions from BBT (2024).

    The Stata code uses:
    - areg for OLS with country FE
    - ivreg2 for IV estimation with cluster-robust SEs
    - cc* for country dummies, yy* for time dummies
    - partial(yy* cc*) to partial out FE in IV

    Paper description (p.728):
    "The first-and second-moment series are scaled for comparability across 
    columns to have residualized unit standard deviation over the regression sample."
    
    This means variables should be scaled so that after partialling out FE,
    they have unit standard deviation.
    """

    def __init__(self, data_path: str = None, standardize_residualized: str = 'none'):
        """
        Initialize Panel IV regressions.
        
        Parameters
        ----------
        data_path : str, optional
            Path to the Stata data file.
        standardize_residualized : str, optional
            How to scale variables:
            - 'none': No scaling (default)
            - 'semi': Scale X to residualized unit std, y unchanged
            - 'full': Scale both X and y to residualized unit std (beta coefficients)
        
        Paper (p.728): "The first-and second-moment series are scaled for 
        comparability across columns to have residualized unit standard deviation 
        over the regression sample."
        """
        if data_path is None:
            data_path = PROJECT_ROOT / "data" / "IV" / "panel_iv_data.dta"
        self.data_path = Path(data_path)
        self.df = None
        self.output_dir = PROJECT_ROOT / "output" / "tables"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # How to scale variables to residualized unit standard deviation
        # 'none': No scaling
        # 'semi': Scale X only (semi-standardized)
        # 'full': Scale both X and y (beta coefficients)
        self.standardize_residualized = standardize_residualized
        self._resid_stds = {}  # Store residualized stds for reporting

        # IV instrument sets (matching Stata globals)
        self.iv = [
            'l1savgnatshock', 'l1savgpolshock',
            'l1savgrevshock', 'l1savgtershock'
        ]
        self.d_iv = [
            'l1savgd_natshock', 'l1savgd_polshock',
            'l1savgd_revshock', 'l1savgd_tershock'
        ]
        self.t_iv = [
            'l1savgt_natshock', 'l1savgt_polshock',
            'l1savgt_revshock', 'l1savgt_tershock'
        ]

    def load_data(self):
        """Load the panel IV dataset."""
        print(f"Loading data from {self.data_path}")
        self.df = pd.read_stata(self.data_path)
        # Convert to float64 for numerical precision
        for col in self.df.columns:
            if self.df[col].dtype in ['float32', 'float64']:
                self.df[col] = self.df[col].astype(np.float64)
        print(f"  Loaded {len(self.df)} observations, {len(self.df.columns)} columns")
        return self

    def _get_fe_arrays(self, df: pd.DataFrame, sample_mask=None):
        """
        Build country and time FE arrays for a given sample.

        Returns cc_cols, yy_cols as numpy arrays.
        """
        cc_cols, yy_cols = get_cc_yy_cols(df)
        cc_arr = df[cc_cols].values.astype(np.float64)
        yy_arr = df[yy_cols].values.astype(np.float64)
        partial = np.hstack([cc_arr, yy_arr])
        return partial, cc_cols, yy_cols

    def _compute_residualized_std(self, x: np.ndarray, partial: np.ndarray) -> float:
        """
        Compute residualized standard deviation.
        
        Paper (p.728): "scaled to have residualized unit standard deviation"
        This means: std(x - Px @ x) where P is the projection matrix for FE.
        """
        # QR decomposition for numerical stability
        Q, R = np.linalg.qr(partial, mode='reduced')
        diag_R = np.abs(np.diag(R))
        rank = np.sum(diag_R > 1e-10 * diag_R[0])
        Q_rank = Q[:, :rank]
        Px = Q_rank @ Q_rank.T
        
        # Residualize
        x_resid = x - Px @ x
        return np.std(x_resid)

    def _scale_to_residualized_unit_std(
        self, 
        x: np.ndarray, 
        partial: np.ndarray,
        var_name: str = None
    ) -> tuple:
        """
        Scale variable to have residualized unit standard deviation.
        
        Returns: (x_scaled, residualized_std)
        """
        rstd = self._compute_residualized_std(x, partial)
        if var_name:
            self._resid_stds[var_name] = rstd
        if rstd > 1e-10:
            return x / rstd, rstd
        else:
            return x, rstd

    def _prepare_iv_regression(
        self,
        df: pd.DataFrame,
        endog_vars: list,
        instr_vars: list,
        cluster_var: str = 'country',
        sample_filter=None,
    ):
        """
        Prepare data for IV regression.

        Handles missing values, builds arrays for iv2sls.
        
        If self.standardize_residualized is True, scales endogenous variables
        to have residualized unit standard deviation (matching paper description).
        """
        data = df.copy()
        if sample_filter is not None:
            data = data[sample_filter].copy()

        # All needed columns
        all_vars = ['ydgdp'] + endog_vars + instr_vars + [cluster_var]
        cc_cols, yy_cols = get_cc_yy_cols(data)
        all_vars += cc_cols + yy_cols

        # Drop rows with missing values in any needed variable
        valid = data[all_vars].notna().all(axis=1)
        data = data[valid].copy()

        y = data['ydgdp'].values.astype(np.float64)
        X_endog = np.column_stack([data[v].values.astype(np.float64) for v in endog_vars])
        Z = np.column_stack([data[v].values.astype(np.float64) for v in instr_vars])
        clusters = data[cluster_var].values

        # Partial out FE (cc* and yy*)
        partial, _, _ = self._get_fe_arrays(data)
        
        # Scale to residualized unit std if requested
        # Paper (p.728): "scaled to have residualized unit standard deviation"
        if self.standardize_residualized in ['semi', 'full']:
            # Scale y (only for 'full' mode)
            if self.standardize_residualized == 'full':
                y, y_rstd = self._scale_to_residualized_unit_std(y, partial, 'ydgdp')
            
            # Scale each endogenous variable (for both 'semi' and 'full')
            X_endog_scaled = np.zeros_like(X_endog)
            for i, var in enumerate(endog_vars):
                X_endog_scaled[:, i], _ = self._scale_to_residualized_unit_std(
                    X_endog[:, i], partial, var
                )
            X_endog = X_endog_scaled

        return y, X_endog, Z, clusters, partial, data

    def _run_iv(
        self,
        endog_vars: list,
        instr_vars: list,
        cluster: bool = True,
        sample_filter=None,
    ):
        """
        Run IV regression matching Stata's ivreg2.

        ivreg2 ydgdp cc* yy* (endog = instr), cluster(country) partial(yy* cc*)
        """
        y, X_endog, Z, clusters, partial, data = self._prepare_iv_regression(
            self.df, endog_vars, instr_vars, sample_filter=sample_filter
        )

        if not cluster:
            # Generate unique cluster per observation for heteroskedastic SEs
            clusters = np.arange(len(y))

        result = iv2sls_with_cluster_se(
            y=y,
            X_endog=X_endog,
            X_exog=np.empty((len(y), 0)),  # No additional exog vars beyond FE
            Z=Z,
            clusters=clusters,
            partial_out=partial,
        )

        return result, data

    def _run_areg(self, sample_filter=None):
        """
        Run OLS with country and period FE matching Stata's areg.

        Stata: areg ydgdp cs_index_ret cs_index_vol i.yq_int, ab(country) cluster(country)

        This includes BOTH:
        - Period FE (i.yq_int) - explicit dummies for each period
        - Country FE (ab(country)) - absorbed (demeaned)

        We use iterative demeaning (Frisch-Waugh-Lovell) for two-way FE.
        """
        data = self.df.copy()
        if sample_filter is not None:
            data = data[sample_filter].copy()

        # Need yq_int for period FE
        vars_needed = ['ydgdp', 'cs_index_ret', 'cs_index_vol', 'country', 'yq_int']
        valid = data[vars_needed].notna().all(axis=1)
        data = data[valid].copy()

        y = data['ydgdp'].values.astype(np.float64)
        X = np.column_stack([
            data['cs_index_ret'].values.astype(np.float64),
            data['cs_index_vol'].values.astype(np.float64),
        ])
        clusters = data['country'].values
        country_groups = data['country'].values
        period_groups = data['yq_int'].values

        # Two-way FE: iterative demeaning (Frisch-Waugh-Lovell)
        # Demean by country first, then by period, iterate until convergence
        unique_countries = np.unique(country_groups)
        unique_periods = np.unique(period_groups)

        max_iter = 100
        tol = 1e-10
        for iteration in range(max_iter):
            y_old = y.copy()
            X_old = X.copy()

            # Demean by country
            for c in unique_countries:
                mask = country_groups == c
                y[mask] -= y[mask].mean()
                X[mask] -= X[mask].mean(axis=0)

            # Demean by period
            for p in unique_periods:
                mask = period_groups == p
                y[mask] -= y[mask].mean()
                X[mask] -= X[mask].mean(axis=0)

            # Check convergence
            if np.max(np.abs(y - y_old)) < tol and np.max(np.abs(X - X_old)) < tol:
                break

        result = ols_with_cluster_se(y, X, clusters)
        return result, data

    # =========================================================================
    # TABLE 1: Descriptive Statistics
    # =========================================================================
    def table1_dstats(self):
        """
        TABLE 1: Descriptive statistics for sample with non-missing ydgdp.

        Matches: estpost summarize ... if ydgdp~=., de
        esttab ... cell("count mean p50 sd min max")
        """
        print("\n" + "=" * 70)
        print("TABLE 1: Descriptive Statistics")
        print("=" * 70)

        vars_desc = [
            'ydgdp', 'cs_index_ret', 'cs_index_vol',
            'avgret', 'lavgvol', 'avgcs_ret', 'lavgcs_vol',
            'ynatshock', 'ysavgnatshock', 'ypolshock', 'ysavgpolshock',
            'yrevshock', 'ysavgrevshock', 'ytershock', 'ysavgtershock',
            'GDP'
        ]

        data = self.df[self.df['ydgdp'].notna()].copy()

        stats = []
        for var in vars_desc:
            if var in data.columns:
                s = data[var]
                stats.append({
                    'Variable': var,
                    'N': int(s.count()),
                    'Mean': s.mean(),
                    'p50': s.median(),
                    'Std': s.std(),
                    'Min': s.min(),
                    'Max': s.max(),
                })

        result_df = pd.DataFrame(stats)
        print(result_df.to_string(index=False, float_format='%.4f'))

        # Save
        out_path = self.output_dir / "table1_dstats.csv"
        result_df.to_csv(out_path, index=False, float_format='%.6f')
        print(f"\n  Saved to {out_path}")

        return result_df

    # =========================================================================
    # TABLE 2: Baseline Results
    # =========================================================================
    def table2_baseline(self):
        """
        TABLE 2: Baseline IV results.

        Columns:
        1. OLS with country FE and period FE (areg with i.yq_int)
        2. IV - micro+macro (cs_index_ret, cs_index_vol)
        3. IV - stock index (l1avgret, l1lavgvol)
        4. IV - stock index, common sample
        5. IV - cross-section (l1avgcs_ret, l1lavgcs_vol)
        """
        print("\n" + "=" * 70)
        print("TABLE 2: Baseline Results")
        print("=" * 70)

        results = {}

        # Col 1: OLS with country FE and period FE
        print("\n--- Column 1: OLS (areg) ---")
        res1, data1 = self._run_areg()
        results['col1_ols'] = res1
        var_names = ['cs_index_ret', 'cs_index_vol']
        print(format_coef_table(
            res1['coef'], res1['se'], var_names,
            title="Col 1: OLS with Country FE and Period FE",
            nobs=res1['nobs'], nclusters=res1['nclusters'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES'}
        ))

        # Col 2: IV - micro+macro
        print("\n--- Column 2: IV (micro+macro) ---")
        res2, data2 = self._run_iv(
            endog_vars=['cs_index_ret', 'cs_index_vol'],
            instr_vars=self.iv,
            cluster=True,
        )
        results['col2_iv'] = res2
        print(format_coef_table(
            res2['coef_endog'], res2['se_endog'],
            ['cs_index_ret', 'cs_index_vol'],
            title="Col 2: IV (micro+macro)",
            nobs=res2['nobs'], nclusters=res2['nclusters'],
            addtext={
                'Period FE': 'YES', 'Country FE': 'YES',
                'Hansen J pval': f"{res2['J_pval']:.4f}" if not np.isnan(res2['J_pval']) else 'N/A',
                'F-Stat (1st)': f"{res2['first_stage'][0]['F_stat']:.4f}",
            }
        ))

        # Col 3: IV - stock index
        print("\n--- Column 3: IV (stock index) ---")
        res3, data3 = self._run_iv(
            endog_vars=['l1avgret', 'l1lavgvol'],
            instr_vars=self.iv,
            cluster=True,
        )
        results['col3_iv_stock'] = res3
        print(format_coef_table(
            res3['coef_endog'], res3['se_endog'],
            ['l1avgret', 'l1lavgvol'],
            title="Col 3: IV (stock index)",
            nobs=res3['nobs'], nclusters=res3['nclusters'],
            addtext={
                'Period FE': 'YES', 'Country FE': 'YES',
                'Hansen J pval': f"{res3['J_pval']:.4f}" if not np.isnan(res3['J_pval']) else 'N/A',
                'F-Stat (1st)': f"{res3['first_stage'][0]['F_stat']:.4f}",
            }
        ))

        # Col 4: IV - stock index, common sample
        print("\n--- Column 4: IV (stock index, common sample) ---")
        # Common sample = sample from Col 1 OLS
        sample1_vars = ['ydgdp', 'cs_index_ret', 'cs_index_vol', 'country']
        common_mask = self.df[sample1_vars].notna().all(axis=1)
        res4, data4 = self._run_iv(
            endog_vars=['l1avgret', 'l1lavgvol'],
            instr_vars=self.iv,
            cluster=True,
            sample_filter=common_mask,
        )
        results['col4_iv_stock_common'] = res4
        print(format_coef_table(
            res4['coef_endog'], res4['se_endog'],
            ['l1avgret', 'l1lavgvol'],
            title="Col 4: IV (stock index, common sample)",
            nobs=res4['nobs'], nclusters=res4['nclusters'],
            addtext={
                'Period FE': 'YES', 'Country FE': 'YES',
                'Hansen J pval': f"{res4['J_pval']:.4f}" if not np.isnan(res4['J_pval']) else 'N/A',
                'F-Stat (1st)': f"{res4['first_stage'][0]['F_stat']:.4f}",
            }
        ))

        # Col 5: IV - cross-section
        print("\n--- Column 5: IV (cross-section) ---")
        res5, data5 = self._run_iv(
            endog_vars=['l1avgcs_ret', 'l1lavgcs_vol'],
            instr_vars=self.iv,
            cluster=True,
        )
        results['col5_iv_cs'] = res5
        print(format_coef_table(
            res5['coef_endog'], res5['se_endog'],
            ['l1avgcs_ret', 'l1lavgcs_vol'],
            title="Col 5: IV (cross-section)",
            nobs=res5['nobs'], nclusters=res5['nclusters'],
            addtext={
                'Period FE': 'YES', 'Country FE': 'YES',
                'Hansen J pval': f"{res5['J_pval']:.4f}" if not np.isnan(res5['J_pval']) else 'N/A',
                'F-Stat (1st)': f"{res5['first_stage'][0]['F_stat']:.4f}",
            }
        ))

        # Save combined results
        self._save_table_results(results, 'table2_baseline')
        return results

    # =========================================================================
    # TABLE 3: Robustness & Higher Moments
    # =========================================================================
    def table3_robustness(self):
        """
        TABLE 3: Robustness checks and higher moments.

        Columns:
        1. Baseline IV
        2. Population weighted
        3. Add skewness
        4. Keep vol and skewness only
        5. HAR adjusted (daily)
        6. HAR adjusted (quarterly)
        """
        print("\n" + "=" * 70)
        print("TABLE 3: Robustness & Higher Moments")
        print("=" * 70)

        results = {}

        # Col 1: Baseline
        print("\n--- Column 1: Baseline ---")
        res, _ = self._run_iv(
            endog_vars=['cs_index_ret', 'cs_index_vol'],
            instr_vars=self.iv,
            cluster=True,  # Stata: cluster(country)
        )
        results['col1_baseline'] = res
        print(format_coef_table(
            res['coef_endog'], res['se_endog'],
            ['cs_index_ret', 'cs_index_vol'],
            title="Col 1: Baseline",
            nobs=res['nobs'], nclusters=res['nclusters'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res['J_pval']:.4f}" if not np.isnan(res['J_pval']) else 'N/A'}
        ))

        # Col 2: Population weighted (analytic weights)
        print("\n--- Column 2: Population Weighted ---")
        res2, _ = self._run_iv_weighted(
            endog_vars=['cs_index_ret', 'cs_index_vol'],
            instr_vars=self.iv,
            weight_var='lpop',
        )
        results['col2_pop_weighted'] = res2
        print(format_coef_table(
            res2['coef_endog'], res2['se_endog'],
            ['cs_index_ret', 'cs_index_vol'],
            title="Col 2: Population Weighted",
            nobs=res2['nobs'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res2['J_pval']:.4f}" if not np.isnan(res2['J_pval']) else 'N/A'}
        ))

        # Col 3: Add skewness
        print("\n--- Column 3: Add Skewness ---")
        res3, _ = self._run_iv(
            endog_vars=['cs_index_ret', 'cs_index_vol', 'cs_index_skew'],
            instr_vars=self.iv,
            cluster=False,
        )
        results['col3_skewness'] = res3
        print(format_coef_table(
            res3['coef_endog'], res3['se_endog'],
            ['cs_index_ret', 'cs_index_vol', 'cs_index_skew'],
            title="Col 3: Add Skewness",
            nobs=res3['nobs'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res3['J_pval']:.4f}" if not np.isnan(res3['J_pval']) else 'N/A'}
        ))

        # Col 4: Keep vol and skewness only
        print("\n--- Column 4: Vol + Skewness Only ---")
        res4, _ = self._run_iv(
            endog_vars=['cs_index_vol', 'cs_index_skew'],
            instr_vars=self.iv,
            cluster=False,
        )
        results['col4_vol_skew'] = res4
        print(format_coef_table(
            res4['coef_endog'], res4['se_endog'],
            ['cs_index_vol', 'cs_index_skew'],
            title="Col 4: Vol + Skewness Only",
            nobs=res4['nobs'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res4['J_pval']:.4f}" if not np.isnan(res4['J_pval']) else 'N/A'}
        ))

        # Col 5: HAR adjusted (daily)
        print("\n--- Column 5: HAR Adjusted (Daily) ---")
        res5, _ = self._run_iv(
            endog_vars=['cs_index_ret_har', 'cs_index_vol_har'],
            instr_vars=self.iv,
            cluster=True,
        )
        results['col5_har_daily'] = res5
        print(format_coef_table(
            res5['coef_endog'], res5['se_endog'],
            ['cs_index_ret_har', 'cs_index_vol_har'],
            title="Col 5: HAR Adjusted (Daily)",
            nobs=res5['nobs'], nclusters=res5['nclusters'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res5['J_pval']:.4f}" if not np.isnan(res5['J_pval']) else 'N/A'}
        ))

        # Col 6: HAR adjusted (quarterly)
        print("\n--- Column 6: HAR Adjusted (Quarterly) ---")
        res6, _ = self._run_iv(
            endog_vars=['cs_index_ret_har_q', 'cs_index_vol_har_q'],
            instr_vars=self.iv,
            cluster=True,
        )
        results['col6_har_quarterly'] = res6
        print(format_coef_table(
            res6['coef_endog'], res6['se_endog'],
            ['cs_index_ret_har_q', 'cs_index_vol_har_q'],
            title="Col 6: HAR Adjusted (Quarterly)",
            nobs=res6['nobs'], nclusters=res6['nclusters'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res6['J_pval']:.4f}" if not np.isnan(res6['J_pval']) else 'N/A'}
        ))

        self._save_table_results(results, 'table3_robustness')
        return results

    def _run_iv_weighted(
        self,
        endog_vars: list,
        instr_vars: list,
        weight_var: str = 'lpop',
        sample_filter=None,
    ):
        """
        IV regression with analytic weights.

        Matches Stata: ivreg2 ... [aw=lpop], cluster(country) partial(yy* cc*)
        """
        data = self.df.copy()
        if sample_filter is not None:
            data = data[sample_filter].copy()

        all_vars = (['ydgdp'] + endog_vars + instr_vars +
                    [weight_var, 'country'])
        cc_cols, yy_cols = get_cc_yy_cols(data)
        all_vars += cc_cols + yy_cols

        valid = data[all_vars].notna().all(axis=1)
        data = data[valid].copy()

        y = data['ydgdp'].values.astype(np.float64)
        X_endog = np.column_stack([data[v].values.astype(np.float64) for v in endog_vars])
        Z = np.column_stack([data[v].values.astype(np.float64) for v in instr_vars])
        w = np.sqrt(data[weight_var].values.astype(np.float64))  # analytic weights
        clusters = data['country'].values

        partial = np.hstack([data[cc_cols].values.astype(np.float64),
                             data[yy_cols].values.astype(np.float64)])

        # Apply weights
        y_w = y * w
        X_endog_w = X_endog * w[:, np.newaxis]
        Z_w = Z * w[:, np.newaxis]
        partial_w = partial * w[:, np.newaxis]

        result = iv2sls_with_cluster_se(
            y=y_w,
            X_endog=X_endog_w,
            X_exog=np.empty((len(y_w), 0)),
            Z=Z_w,
            clusters=clusters,
            partial_out=partial_w,
        )

        return result, data

    # =========================================================================
    # TABLE 4: Trade and Distance Weighted Instruments
    # =========================================================================
    def table4_weighting(self):
        """
        TABLE 4: Trade and distance weighted instruments.

        Uses common sample from Table 2 Col 1.
        """
        print("\n" + "=" * 70)
        print("TABLE 4: Trade and Distance Weighted Instruments")
        print("=" * 70)

        results = {}

        # Common sample
        sample1_vars = ['ydgdp', 'cs_index_ret', 'cs_index_vol', 'country']
        common_mask = self.df[sample1_vars].notna().all(axis=1)

        specs = [
            ('col1_trade_pcf', ['cs_index_ret', 'cs_index_vol'], self.t_iv),
            ('col2_dist_pcf', ['cs_index_ret', 'cs_index_vol'], self.d_iv),
            ('col3_trade_stock', ['l1avgret', 'l1lavgvol'], self.t_iv),
            ('col4_dist_stock', ['l1avgret', 'l1lavgvol'], self.d_iv),
            ('col5_trade_cs', ['l1avgcs_ret', 'l1lavgcs_vol'], self.t_iv),
            ('col6_dist_cs', ['l1avgcs_ret', 'l1lavgcs_vol'], self.d_iv),
        ]

        for name, endog, instr in specs:
            print(f"\n--- {name} ---")
            res, _ = self._run_iv(
                endog_vars=endog,
                instr_vars=instr,
                cluster=True,
                sample_filter=common_mask,
            )
            results[name] = res
            print(format_coef_table(
                res['coef_endog'], res['se_endog'], endog,
                title=name,
                nobs=res['nobs'], nclusters=res['nclusters'],
                addtext={'Period FE': 'YES', 'Country FE': 'YES',
                         'Hansen J pval': f"{res['J_pval']:.4f}" if not np.isnan(res['J_pval']) else 'N/A'}
            ))

        self._save_table_results(results, 'table4_weighting')
        return results

    # =========================================================================
    # TABLE 5: Media Weightings and Jump Definitions
    # =========================================================================
    def table5_media_weightings(self):
        """
        TABLE 5: Alternative media weightings and jump definitions.
        """
        print("\n" + "=" * 70)
        print("TABLE 5: Media Weightings and Jump Definitions")
        print("=" * 70)

        results = {}

        specs = [
            ('col1_baseline', ['cs_index_ret', 'cs_index_vol'], self.iv),
            ('col2_unscaled', ['cs_index_ret', 'cs_index_vol'],
             ['l1s0avgnatshock', 'l1s0avgpolshock', 'l1s0avgrevshock', 'l1s0avgtershock']),
            ('col3_over_median', ['cs_index_ret', 'cs_index_vol'],
             ['l1sMedavgnatshock', 'l1sMedavgpolshock', 'l1sMedavgrevshock', 'l1sMedavgtershock']),
            ('col4_narrower_window', ['cs_index_ret', 'cs_index_vol'],
             ['l1s2avgnatshock', 'l1s2avgpolshock', 'l1s2avgrevshock', 'l1s2avgtershock']),
            ('col5_scaled_jump_reg', ['cs_index_ret', 'cs_index_vol'],
             ['l1sprdavgnatshock', 'l1sprdavgpolshock', 'l1sprdavgrevshock', 'l1sprdavgtershock']),
            ('col6_nonwestern', ['cs_index_ret', 'cs_index_vol'],
             ['l1nwsavgnatshock', 'l1nwsavgpolshock', 'l1nwsavgrevshock', 'l1nwsavgtershock']),
        ]

        for name, endog, instr in specs:
            print(f"\n--- {name} ---")
            res, _ = self._run_iv(
                endog_vars=endog,
                instr_vars=instr,
                cluster=True,
            )
            results[name] = res
            print(format_coef_table(
                res['coef_endog'], res['se_endog'], endog,
                title=name,
                nobs=res['nobs'], nclusters=res['nclusters'],
                addtext={'Period FE': 'YES', 'Country FE': 'YES',
                         'Hansen J pval': f"{res['J_pval']:.4f}" if not np.isnan(res['J_pval']) else 'N/A'}
            ))

        self._save_table_results(results, 'table5_media_weightings')
        return results

    # =========================================================================
    # TABLE 6: Alternative Uncertainty Proxies
    # =========================================================================
    def table6_alternative_uncertainty(self):
        """
        TABLE 6: Alternative uncertainty proxies.

        Columns:
        1. WUI
        2. EPU
        3. EPU + WUI
        4. Consensus Forecasts
        5. Exchange Rates
        """
        print("\n" + "=" * 70)
        print("TABLE 6: Alternative Uncertainty Proxies")
        print("=" * 70)

        results = {}

        # Col 1: WUI
        print("\n--- Column 1: WUI ---")
        res1, _ = self._run_iv(
            endog_vars=['cs_index_ret', 'l1lavgWUI'],
            instr_vars=self.iv,
            cluster=False,
        )
        results['col1_wui'] = res1
        print(format_coef_table(
            res1['coef_endog'], res1['se_endog'],
            ['cs_index_ret', 'l1lavgWUI'],
            title="Col 1: WUI",
            nobs=res1['nobs'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res1['J_pval']:.4f}" if not np.isnan(res1['J_pval']) else 'N/A'}
        ))

        # Col 2: EPU
        print("\n--- Column 2: EPU ---")
        res2, _ = self._run_iv(
            endog_vars=['cs_index_ret', 'l1lavgEPU'],
            instr_vars=self.iv,
            cluster=False,
        )
        results['col2_epu'] = res2
        print(format_coef_table(
            res2['coef_endog'], res2['se_endog'],
            ['cs_index_ret', 'l1lavgEPU'],
            title="Col 2: EPU",
            nobs=res2['nobs'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res2['J_pval']:.4f}" if not np.isnan(res2['J_pval']) else 'N/A'}
        ))

        # Col 3: EPU + WUI (combined)
        print("\n--- Column 3: EPU + WUI ---")
        # Create combined variable: use EPU where available, else WUI
        data = self.df.copy()
        any_epu = data.groupby('country')['EPU'].transform('mean')
        epu_wui = data['l1lavgEPU'].copy()
        wui_mask = data['l1lavgEPU'].isna() & any_epu.isna()
        epu_wui[wui_mask] = data.loc[wui_mask, 'l1lavgWUI']
        data['epu_wui'] = epu_wui

        # Need to manually run IV with this derived variable
        all_vars = (['ydgdp', 'cs_index_ret', 'epu_wui'] + self.iv +
                    ['country'])
        cc_cols, yy_cols = get_cc_yy_cols(data)
        all_vars += cc_cols + yy_cols
        valid = data[all_vars].notna().all(axis=1)
        data_valid = data[valid].copy()

        y = data_valid['ydgdp'].values.astype(np.float64)
        X_endog = np.column_stack([
            data_valid['cs_index_ret'].values.astype(np.float64),
            data_valid['epu_wui'].values.astype(np.float64),
        ])
        Z = np.column_stack([data_valid[v].values.astype(np.float64) for v in self.iv])
        clusters = data_valid['country'].values
        partial = np.hstack([data_valid[cc_cols].values.astype(np.float64),
                             data_valid[yy_cols].values.astype(np.float64)])

        res3 = iv2sls_with_cluster_se(
            y=y, X_endog=X_endog, X_exog=np.empty((len(y), 0)),
            Z=Z, clusters=clusters, partial_out=partial,
        )
        results['col3_epu_wui'] = res3
        print(format_coef_table(
            res3['coef_endog'], res3['se_endog'],
            ['cs_index_ret', 'epu_wui'],
            title="Col 3: EPU + WUI",
            nobs=res3['nobs'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res3['J_pval']:.4f}" if not np.isnan(res3['J_pval']) else 'N/A'}
        ))

        # Col 4: Consensus Forecasts
        print("\n--- Column 4: Consensus Forecasts ---")
        res4, _ = self._run_iv(
            endog_vars=['cs_index_ret', 'l1lgdp_for_sd'],
            instr_vars=self.iv,
            cluster=False,
        )
        results['col4_consensus'] = res4
        print(format_coef_table(
            res4['coef_endog'], res4['se_endog'],
            ['cs_index_ret', 'l1lgdp_for_sd'],
            title="Col 4: Consensus Forecasts",
            nobs=res4['nobs'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res4['J_pval']:.4f}" if not np.isnan(res4['J_pval']) else 'N/A'}
        ))

        # Col 5: Exchange Rates
        print("\n--- Column 5: Exchange Rates ---")
        res5, _ = self._run_iv(
            endog_vars=['cs_index_ret', 'l1lavgexchgvol'],
            instr_vars=self.iv,
            cluster=True,
        )
        results['col5_exchg'] = res5
        print(format_coef_table(
            res5['coef_endog'], res5['se_endog'],
            ['cs_index_ret', 'l1lavgexchgvol'],
            title="Col 5: Exchange Rates",
            nobs=res5['nobs'], nclusters=res5['nclusters'],
            addtext={'Period FE': 'YES', 'Country FE': 'YES',
                     'Hansen J pval': f"{res5['J_pval']:.4f}" if not np.isnan(res5['J_pval']) else 'N/A'}
        ))

        self._save_table_results(results, 'table6_alternative_uncertainty')
        return results

    # =========================================================================
    # Utility
    # =========================================================================
    def _save_table_results(self, results: dict, table_name: str):
        """Save regression results to CSV."""
        rows = []
        for col_name, res in results.items():
            if 'coef_endog' in res:
                # IV results
                n_endog = len(res['coef_endog'])
                for i in range(n_endog):
                    row = {
                        'specification': col_name,
                        'variable': f'endog_{i}',
                        'coefficient': res['coef_endog'][i],
                        'std_error': res['se_endog'][i],
                        'nobs': res['nobs'],
                        'J_stat': res['J_stat'],
                        'J_pval': res['J_pval'],
                    }
                    if i < len(res['first_stage']):
                        row['F_stat_1st'] = res['first_stage'][i]['F_stat']
                    rows.append(row)
            elif 'coef' in res:
                # OLS results
                for i in range(len(res['coef'])):
                    rows.append({
                        'specification': col_name,
                        'variable': f'exog_{i}',
                        'coefficient': res['coef'][i],
                        'std_error': res['se'][i],
                        'nobs': res['nobs'],
                        'J_stat': np.nan,
                        'J_pval': np.nan,
                    })

        df = pd.DataFrame(rows)
        out_path = self.output_dir / f"{table_name}.csv"
        df.to_csv(out_path, index=False, float_format='%.6f')
        print(f"\n  Saved to {out_path}")

    def run_all(self):
        """Run all tables."""
        self.load_data()
        self.table1_dstats()
        self.table2_baseline()
        self.table3_robustness()
        self.table4_weighting()
        self.table5_media_weightings()
        self.table6_alternative_uncertainty()
        print("\n" + "=" * 70)
        print("ALL TABLES COMPLETE")
        print("=" * 70)


if __name__ == '__main__':
    panel_iv = PanelIV()
    panel_iv.run_all()
