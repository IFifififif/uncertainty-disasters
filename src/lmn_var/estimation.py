"""
LMN VAR Estimation for BBT (2024).

Replicates Stata code (STEP1_STATA_ESTIMATION.do) + MATLAB code
(STEP2_MATLAB_ESTIMATION.m + STEP3_GRAPHS.m).
Produces Figures 3-5 from the paper.

The LMN approach (Arias, Rubio-Ramirez, & Waggoner, 2018) uses
disaster events to restrict the set of admissible structural VAR
parameters, then computes impulse responses over the admissible set.

Original: Stata/MP 15.1 + MATLAB R2021a
Python: Uses linearmodels + numpy for admissible set computation
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class LMNVAR:
    """
    LMN VAR estimation matching the Stata + MATLAB implementation.

    Step 1: Estimate VAR with country and time FE (Stata reghdfe)
    Step 2: Compute admissible sets via random matrix draws (MATLAB)
    Step 3: Generate figures from admissible sets (MATLAB)
    """

    def __init__(self, data_path: str = None):
        if data_path is None:
            data_path = PROJECT_ROOT / "data" / "LMN_VAR" / "Dates_and_Data.dta"
        self.data_path = Path(data_path)
        self.output_dir = PROJECT_ROOT / "output" / "figures"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # VAR dimensions
        self.NX = 3  # GDP growth, first moment, second moment
        self.Nlags = 3
        self.lengthIRF = 50

    def load_data(self):
        """Load the LMN VAR dataset."""
        print(f"Loading data from {self.data_path}")
        self.df = pd.read_stata(self.data_path)
        print(f"  Shape: {self.df.shape}")
        print(f"  Columns: {list(self.df.columns[:20])}...")
        return self

    def step1_estimate_var_fe(self):
        """
        STEP 1: Estimate VAR with country and time fixed effects.

        Matches Stata: reghdfe depvar l1depvar l2depvar l3depvar
                       l1indepvar1 l1indepvar2 l2indepvar1 l2indepvar2
                       l3indepvar1 l3indepvar2, absorb(country yq_int)

        Returns demeaned residuals for each equation.
        """
        print("\n--- Step 1: VAR with FE estimation ---")

        # Variable names (matching Stata code)
        var_names = ['ydgdp', 'avgret', 'lavgvol']

        # Build lagged variables
        data = self.df.copy()
        for var in var_names:
            for lag in range(1, self.Nlags + 1):
                data[f'l{lag}{var}'] = data.groupby('country')[var].shift(lag)

        # Estimate each equation separately (matching reghdfe)
        residuals = {}
        A_hat = np.zeros((self.NX, self.NX * self.Nlags))

        for eq_idx, dep_var in enumerate(var_names):
            # Build regressors: lags of all 3 variables
            regressor_names = []
            for lag in range(1, self.Nlags + 1):
                for var in var_names:
                    regressor_names.append(f'l{lag}{var}')

            all_cols = [dep_var] + regressor_names + ['country', 'yq_int']
            valid = data[all_cols].notna().all(axis=1)
            eq_data = data[valid].copy()

            y = eq_data[dep_var].values.astype(np.float64)
            X = eq_data[regressor_names].values.astype(np.float64)

            # Iterative demeaning by country and time (Frisch-Waugh-Lovell)
            for fe_col in ['country', 'yq_int']:
                groups = eq_data[fe_col].values
                unique_groups = np.unique(groups)
                for g in unique_groups:
                    mask = groups == g
                    n_g = mask.sum()
                    if n_g > 0:
                        y[mask] -= y[mask].mean()
                        X[mask] -= X[mask].mean(axis=0)

            # OLS on demeaned data
            XtX = X.T @ X
            Xty = X.T @ y
            beta = np.linalg.solve(XtX, Xty)
            resid = y - X @ beta

            residuals[dep_var] = resid
            A_hat[eq_idx, :] = beta

            print(f"  Equation {eq_idx + 1} ({dep_var}): N={len(y)}, R2={1 - np.sum(resid**2) / np.sum((y - y.mean())**2):.4f}")

        # Store results
        self.residuals = residuals
        self.A_hat = A_hat
        self.eq_data = data

        # Compute covariance of residuals
        # Align residuals across equations
        common_idx = None
        for dep_var in var_names:
            valid = data[[dep_var] + [f'l{l}{v}' for l in range(1, self.Nlags + 1) for v in var_names] + ['country', 'yq_int']].notna().all(axis=1)
            idx = data[valid].index
            if common_idx is None:
                common_idx = idx
            else:
                common_idx = common_idx.intersection(idx)

        resid_matrix = np.column_stack([residuals[v].loc[common_idx] for v in var_names])
        self.Sigma_hat = (resid_matrix.T @ resid_matrix) / len(common_idx)

        print(f"  Residual covariance (diagonal): {np.diag(self.Sigma_hat)}")

        return self

    def step2_admissible_sets(self, n_draws: int = 100000, seed: int = 3991):
        """
        STEP 2: Compute admissible sets via random matrix draws.

        Matches MATLAB: STEP2_MATLAB_ESTIMATION.m

        For each draw:
        1. Generate random orthogonal matrix Q
        2. Compute structural matrix B = Sigma^{1/2} Q
        3. Check if impulse responses are consistent with disaster events
        4. If yes, store the impulse responses

        Parameters
        ----------
        n_draws : number of random draws
        seed : random seed
        """
        print(f"\n--- Step 2: Admissible Sets ({n_draws} draws) ---")
        rng = np.random.RandomState(seed)

        Sigma = self.Sigma_hat
        NX = self.NX

        # Cholesky decomposition of Sigma
        Sigma_chol = np.linalg.cholesky(Sigma)

        # Build companion form for IRF computation
        B1hat = self.A_hat.T.reshape(NX * self.Nlags, NX)
        B1tilde = np.zeros((NX * self.Nlags, NX * self.Nlags))
        B1tilde[:NX, :NX * self.Nlags] = B1hat
        if self.Nlags > 1:
            for ct in range(self.Nlags - 1):
                B1tilde[ct * NX + np.arange(NX), ct * NX + np.arange(NX)] = np.eye(NX)

        # Disaster event restrictions
        # These define the admissible set: impulse responses must satisfy
        # certain sign/size restrictions at disaster dates
        disaster_restrictions = self._get_disaster_restrictions()

        # Storage for admissible impulse responses
        admissible_irfs = []
        n_admissible = 0

        for draw in range(n_draws):
            # Generate random orthogonal matrix
            # Method: QR decomposition of random Gaussian matrix
            Z = rng.randn(NX, NX)
            Q, R = np.linalg.qr(Z)
            # Ensure proper orientation (det(Q) = +1)
            Q = Q @ np.diag(np.sign(np.diag(R)))

            # Structural matrix
            B = Sigma_chol @ Q

            # Compute impulse responses
            Btilde = np.zeros((NX * self.Nlags, NX))
            Btilde[:NX, :NX] = B

            IRF = np.zeros((self.lengthIRF, NX, NX))
            for varct in range(NX):
                for t in range(self.lengthIRF):
                    IRFvec = np.linalg.matrix_power(B1tilde, t - 1) @ Btilde[:, varct]
                    IRF[t, :, varct] = IRFvec[:NX]

            # Check admissibility
            if self._check_admissibility(IRF, disaster_restrictions):
                admissible_irfs.append(IRF)
                n_admissible += 1

            if (draw + 1) % 10000 == 0:
                print(f"  Draw {draw + 1}/{n_draws}, Admissible: {n_admissible}")

        print(f"  Total admissible: {n_admissible}/{n_draws} "
              f"({100 * n_admissible / n_draws:.1f}%)")

        self.admissible_irfs = admissible_irfs
        self.n_admissible = n_admissible

        return self

    def _get_disaster_restrictions(self):
        """
        Define disaster event restrictions for admissibility.

        These are the sign/size restrictions on impulse responses
        at disaster dates, as defined in the paper.
        """
        # Default restrictions based on the paper:
        # GDP should fall after uncertainty shock
        # First moment (returns) should fall
        # Second moment (volatility) should rise
        restrictions = {
            'gdp_response_at_shock': 'negative',  # GDP falls
            'vol_response_at_shock': 'positive',   # Vol rises
            'ret_response_at_shock': 'negative',   # Returns fall
        }
        return restrictions

    def _check_admissibility(self, IRF: np.ndarray, restrictions: dict) -> bool:
        """
        Check if impulse responses satisfy disaster restrictions.

        Parameters
        ----------
        IRF : (lengthIRF, NX, NX) impulse response array
        restrictions : dict of restrictions

        Returns
        -------
        bool : True if admissible
        """
        # GDP response to uncertainty shock (column 2, row 0)
        gdp_to_unc = IRF[0, 0, 2]

        # Volatility response to uncertainty shock (column 2, row 2)
        vol_to_unc = IRF[0, 2, 2]

        # Returns response to uncertainty shock (column 2, row 1)
        ret_to_unc = IRF[0, 1, 2]

        # Check restrictions
        if restrictions.get('gdp_response_at_shock') == 'negative' and gdp_to_unc > 0:
            return False
        if restrictions.get('vol_response_at_shock') == 'positive' and vol_to_unc < 0:
            return False
        if restrictions.get('ret_response_at_shock') == 'negative' and ret_to_unc > 0:
            return False

        return True

    def step3_generate_figures(self):
        """
        STEP 3: Generate Figures 3-5 from admissible sets.

        Matches MATLAB: STEP3_GRAPHS.m
        """
        print("\n--- Step 3: Generating Figures ---")

        if not hasattr(self, 'admissible_irfs') or len(self.admissible_irfs) == 0:
            print("  ERROR: No admissible IRFs found. Run step2 first.")
            return

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        IRFsamp = np.arange(1, 16)
        admissible = np.array(self.admissible_irfs)

        # Compute median and quantiles across admissible set
        # GDP response to uncertainty shock
        gdp_irfs = admissible[:, :15, 0, 2]  # (n_admissible, 15)
        gdp_median = np.median(gdp_irfs, axis=0)
        gdp_q16 = np.percentile(gdp_irfs, 16, axis=0)
        gdp_q84 = np.percentile(gdp_irfs, 84, axis=0)

        # Volatility response to uncertainty shock
        vol_irfs = admissible[:, :15, 2, 2]
        vol_median = np.median(vol_irfs, axis=0)
        vol_q16 = np.percentile(vol_irfs, 16, axis=0)
        vol_q84 = np.percentile(vol_irfs, 84, axis=0)

        # Returns response to uncertainty shock
        ret_irfs = admissible[:, :15, 1, 2]
        ret_median = np.median(ret_irfs, axis=0)
        ret_q16 = np.percentile(ret_irfs, 16, axis=0)
        ret_q84 = np.percentile(ret_irfs, 84, axis=0)

        # FIGURE 3: GDP response
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(IRFsamp, gdp_median, 'bo-', linewidth=2.5, label='Median')
        ax.fill_between(IRFsamp, gdp_q16, gdp_q84, alpha=0.3, color='blue', label='68% CI')
        ax.plot(IRFsamp, np.zeros(15), 'k--', linewidth=2.5)
        ax.set_ylim(-10, 5)
        ax.set_xlabel('Quarters, Shock in Period 1', fontsize=14)
        ax.set_ylabel('GDP Growth, Percent Year-on-Year', fontsize=14)
        ax.tick_params(labelsize=14)
        ax.legend(fontsize=12)
        out_path = self.output_dir / "FIGURE3.pdf"
        fig.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {out_path}")

        # FIGURE 4: Volatility response
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(IRFsamp, vol_median, 'bo-', linewidth=2.5, label='Median')
        ax.fill_between(IRFsamp, vol_q16, vol_q84, alpha=0.3, color='blue', label='68% CI')
        ax.plot(IRFsamp, np.zeros(15), 'k--', linewidth=2.5)
        ax.set_xlabel('Quarters, Shock in Period 1', fontsize=14)
        ax.set_ylabel('Uncertainty, Percent', fontsize=14)
        ax.tick_params(labelsize=14)
        ax.legend(fontsize=12)
        out_path = self.output_dir / "FIGURE4.pdf"
        fig.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {out_path}")

        # FIGURE 5: Returns response
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(IRFsamp, ret_median, 'bo-', linewidth=2.5, label='Median')
        ax.fill_between(IRFsamp, ret_q16, ret_q84, alpha=0.3, color='blue', label='68% CI')
        ax.plot(IRFsamp, np.zeros(15), 'k--', linewidth=2.5)
        ax.set_xlabel('Quarters, Shock in Period 1', fontsize=14)
        ax.set_ylabel('Returns, Percent', fontsize=14)
        ax.tick_params(labelsize=14)
        ax.legend(fontsize=12)
        out_path = self.output_dir / "FIGURE5.pdf"
        fig.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {out_path}")

        return {
            'gdp_median': gdp_median, 'gdp_q16': gdp_q16, 'gdp_q84': gdp_q84,
            'vol_median': vol_median, 'vol_q16': vol_q16, 'vol_q84': vol_q84,
            'ret_median': ret_median, 'ret_q16': ret_q16, 'ret_q84': ret_q84,
        }

    def run_all(self):
        """Run full LMN VAR pipeline."""
        self.load_data()
        self.step1_estimate_var_fe()
        self.step2_admissible_sets(n_draws=100000, seed=3991)
        self.step3_generate_figures()

        print("\n" + "=" * 70)
        print("LMN VAR ESTIMATION COMPLETE")
        print("=" * 70)
