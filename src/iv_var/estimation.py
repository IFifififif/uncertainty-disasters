"""
IV-VAR Estimation for BBT (2024).

Replicates MATLAB code: STEP1_ESTIMATION.m + STEP2_GRAPHS.m
Produces Figures 6-7 from the paper.

Original: MATLAB R2021a, ~1 min runtime, ~1.5 GB memory
Python: Uses scipy.optimize.minimize for GMM, numpy for IRF computation
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize
from typing import Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class IVVAR:
    """
    IV-VAR estimation matching the MATLAB implementation.

    The VAR has 3 variables (NX=3):
      1. GDP growth (ydgdp)
      2. First moment proxy (avgret)
      3. Second moment proxy (lavgvol)

    With 4 instruments (ND=4):
      1. Natural disaster shock
      2. Political shock
      3. Revolution shock
      4. Terrorist attack shock

    Parameters estimated (Nparams=17):
      - B matrix (3x3 = 9 params): contemporaneous impact matrix
      - Dcoeff matrix (4x2 = 8 params): IV first-stage coefficients
    """

    def __init__(self, data_path: str = None):
        if data_path is None:
            data_path = PROJECT_ROOT / "data" / "IV_VAR" / "VARdata.csv"
        self.data_path = Path(data_path)
        self.output_dir = PROJECT_ROOT / "output" / "figures"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Dimensions
        self.NX = 3  # number of variables in VAR
        self.ND = 4  # number of instruments
        self.Nparams = 17  # 9 (B) + 8 (Dcoeff)
        self.Nmoms = 18  # 6 (Omega) + 12 (E[D*eta])

        # Default settings
        self.Nlags = 3
        self.lengthIRF = 50

    def load_data(self):
        """Load VARdata.csv."""
        print(f"Loading data from {self.data_path}")
        self.data = pd.read_csv(self.data_path)
        print(f"  Shape: {self.data.shape}")
        return self

    def _build_moment_vector(self, X: np.ndarray, D: np.ndarray):
        """
        Build the empirical moment vector from data.

        Parameters
        ----------
        X : (T, NX) VAR data matrix
        D : (T, ND) instrument matrix

        Returns
        -------
        MOMvec : (Nmoms,) empirical moment vector
        """
        T = X.shape[0]

        # Estimate VAR coefficients via OLS
        # X_t = A1 X_{t-1} + ... + A_p X_{t-p} + eta_t
        Y = X[self.Nlags:]  # (T-p, NX)
        n_lag_obs = T - self.Nlags

        # Build lagged regressor matrix
        X_lagged = np.zeros((n_lag_obs, self.NX * self.Nlags))
        for lag in range(1, self.Nlags + 1):
            X_lagged[:, (lag - 1) * self.NX: lag * self.NX] = X[self.Nlags - lag: T - lag]

        # OLS: eta_hat = Y - X_lagged @ A_hat
        # A_hat = (X_lagged' X_lagged)^{-1} X_lagged' Y
        XtX = X_lagged.T @ X_lagged
        XtY = X_lagged.T @ Y
        A_hat = np.linalg.solve(XtX, XtY)
        eta_hat = Y - X_lagged @ A_hat  # (T-p, NX)

        # Covariance of reduced-form residuals: Omega = E[eta eta']
        Omega = (eta_hat.T @ eta_hat) / n_lag_obs

        # E[D_t * eta_t] for each instrument
        D_lagged = D[self.Nlags - 1: T - 1]  # align with eta_hat
        EDeta = np.zeros((self.NX, self.ND))
        for iv_ct in range(self.ND):
            EDeta[:, iv_ct] = (D_lagged[:, iv_ct][:, np.newaxis] * eta_hat).sum(axis=0) / n_lag_obs

        # Build moment vector
        MOMvec = np.zeros(self.Nmoms)
        # Omega moments (6 unique elements of 3x3 symmetric matrix)
        MOMvec[0] = Omega[0, 0]
        MOMvec[1] = Omega[1, 1]
        MOMvec[2] = Omega[2, 2]
        MOMvec[3] = Omega[0, 1]
        MOMvec[4] = Omega[0, 2]
        MOMvec[5] = Omega[1, 2]
        # E[D*eta] moments (3 x 4 = 12)
        for iv_ct in range(self.ND):
            MOMvec[6 + iv_ct * 3 + 0] = EDeta[0, iv_ct]
            MOMvec[6 + iv_ct * 3 + 1] = EDeta[1, iv_ct]
            MOMvec[6 + iv_ct * 3 + 2] = EDeta[2, iv_ct]

        return MOMvec, eta_hat, A_hat

    def _gmm_objective(self, x: np.ndarray, MOMvec: np.ndarray, extraoutput: int = 0):
        """
        GMM objective function (matches MATLAB's fGMMobj.m).

        Parameters
        ----------
        x : (Nparams,) parameter vector
        MOMvec : (Nmoms,) empirical moment vector
        extraoutput : if 1, return implied moments instead of objective

        Returns
        -------
        GMM objective value (sum of squared moment errors)
        """
        NX, ND = self.NX, self.ND

        # Extract B matrix (3x3)
        B = np.zeros((NX, NX))
        B[0, 0] = x[0]; B[1, 0] = x[1]; B[2, 0] = x[2]
        B[0, 1] = x[3]; B[1, 1] = x[4]; B[2, 1] = x[5]
        B[0, 2] = x[6]; B[1, 2] = x[7]; B[2, 2] = x[8]

        # Extract Dcoeff matrix (4x2)
        Dcoeff = np.zeros((ND, 2))
        Dcoeff[0, 0] = x[9];  Dcoeff[1, 0] = x[10]
        Dcoeff[2, 0] = x[11]; Dcoeff[3, 0] = x[12]
        Dcoeff[0, 1] = x[13]; Dcoeff[1, 1] = x[14]
        Dcoeff[2, 1] = x[15]; Dcoeff[3, 1] = x[16]

        # Implied Lambda matrix
        LAMBDA = np.array([
            [1, 0, 0],
            [0, np.sum(Dcoeff[:, 0] ** 2) + 1, np.sum(Dcoeff[:, 0] * Dcoeff[:, 1])],
            [0, np.sum(Dcoeff[:, 0] * Dcoeff[:, 1]), np.sum(Dcoeff[:, 1] ** 2) + 1],
        ])

        # Implied covariance: Omega = B * Lambda * B'
        OMEGA = B @ LAMBDA @ B.T

        # Implied E[D*eta]
        EDeta = np.zeros((NX, ND))
        for iv_ct in range(ND):
            EDeta[:, iv_ct] = B @ np.array([0, Dcoeff[iv_ct, 0], Dcoeff[iv_ct, 1]])

        # Build implied moment vector
        MOMimplied = np.zeros(self.Nmoms)
        MOMimplied[0] = OMEGA[0, 0]
        MOMimplied[1] = OMEGA[1, 1]
        MOMimplied[2] = OMEGA[2, 2]
        MOMimplied[3] = OMEGA[0, 1]
        MOMimplied[4] = OMEGA[0, 2]
        MOMimplied[5] = OMEGA[1, 2]
        for iv_ct in range(ND):
            MOMimplied[6 + iv_ct * 3 + 0] = EDeta[0, iv_ct]
            MOMimplied[6 + iv_ct * 3 + 1] = EDeta[1, iv_ct]
            MOMimplied[6 + iv_ct * 3 + 2] = EDeta[2, iv_ct]

        if extraoutput == 1:
            return MOMimplied

        GMMerr = MOMvec - MOMimplied
        return np.sum(GMMerr ** 2)

    def _compute_irf(self, Bhat: np.ndarray, B1hat: np.ndarray,
                     X: np.ndarray, lengthIRF: int = 50):
        """
        Compute impulse response functions.

        Matches MATLAB IRF computation in VAR.m.
        """
        NX = self.NX
        Nlags = self.Nlags

        # Build companion form
        B1tilde = np.zeros((NX * Nlags, NX * Nlags))
        B1tilde[:NX, :NX * Nlags] = B1hat
        if Nlags > 1:
            for ct in range(Nlags - 1):
                row_start = ct * NX
                col_start = ct * NX
                B1tilde[row_start:row_start + NX, col_start:col_start + NX] = np.eye(NX)

        Btilde = np.zeros((NX * Nlags, NX))
        Btilde[:NX, :NX] = Bhat

        # Compute IRFs
        IRF = np.zeros((lengthIRF, NX, NX))
        for varct in range(NX):
            for t in range(lengthIRF):
                IRFvec = np.linalg.matrix_power(B1tilde, t - 1) @ Btilde[:, varct]
                IRF[t, :, varct] = IRFvec[:NX]
            # Scale: sqrt(var(X(:,varct))) * IRF / Bhat(varct,varct)
            IRF[:, :, varct] = (np.sqrt(np.var(X[:, varct])) *
                                IRF[:, :, varct] / Bhat[varct, varct])

        return IRF

    def _initial_params(self):
        """Return initial parameter guess (matching MATLAB)."""
        p = np.zeros(self.Nparams)
        p[0] = 1; p[5] = 1; p[8] = 1  # B diagonal
        p[9:13] = -1  # Dcoeff col 1
        p[14:17] = 1  # Dcoeff col 2 (first 3)
        p[13] = 0    # Dcoeff(4,1)
        return 0.25 * p

    def estimate_baseline(self, seed: int = 3991):
        """
        Estimate the baseline IV-VAR (matches BASELINE/VAR.m).

        Parameters
        ----------
        seed : random seed (MATLAB uses 3991)

        Returns
        -------
        dict with IRF_S_TO_Y, Bhat, Dcoeffhat, etc.
        """
        print("\n--- Baseline IV-VAR Estimation ---")
        np.random.seed(seed)

        # Load and prepare data
        X = self.data.values[:, :self.NX].astype(np.float64)
        D = self.data.values[:, self.NX:self.NX + self.ND].astype(np.float64)
        T = X.shape[0]

        # Build moments
        MOMvec, eta_hat, A_hat = self._build_moment_vector(X, D)

        # Build VAR coefficient matrix B1hat
        B1hat = A_hat.T.reshape(self.NX, self.NX * self.Nlags)

        # Optimize GMM objective
        param0 = self._initial_params()
        print(f"  Initial GMM objective: {self._gmm_objective(param0, MOMvec):.6f}")

        result = minimize(
            self._gmm_objective,
            param0,
            args=(MOMvec,),
            method='L-BFGS-B',
            options={'maxiter': 50000, 'ftol': 1e-15, 'gtol': 1e-10},
        )

        paramhat = result.x
        print(f"  Final GMM objective: {result.fun:.10f}")
        print(f"  Converged: {result.success}")

        # Extract B matrix
        Bhat = np.zeros((self.NX, self.NX))
        Bhat[0, 0] = paramhat[0]; Bhat[1, 0] = paramhat[1]; Bhat[2, 0] = paramhat[2]
        Bhat[0, 1] = paramhat[3]; Bhat[1, 1] = paramhat[4]; Bhat[2, 1] = paramhat[5]
        Bhat[0, 2] = paramhat[6]; Bhat[1, 2] = paramhat[7]; Bhat[2, 2] = paramhat[8]

        # Extract Dcoeff
        Dcoeffhat = np.zeros((self.ND, 2))
        Dcoeffhat[0, 0] = paramhat[9];  Dcoeffhat[1, 0] = paramhat[10]
        Dcoeffhat[2, 0] = paramhat[11]; Dcoeffhat[3, 0] = paramhat[12]
        Dcoeffhat[0, 1] = paramhat[13]; Dcoeffhat[1, 1] = paramhat[14]
        Dcoeffhat[2, 1] = paramhat[15]; Dcoeffhat[3, 1] = paramhat[16]

        # Compute IRFs
        IRF = self._compute_irf(Bhat, B1hat, X, self.lengthIRF)
        IRF_S_TO_Y = IRF[:15, 0, 2]  # GDP response to uncertainty shock

        print(f"  IRF_S_TO_Y (first 5): {IRF_S_TO_Y[:5]}")

        return {
            'IRF': IRF,
            'IRF_S_TO_Y': IRF_S_TO_Y,
            'Bhat': Bhat,
            'Dcoeffhat': Dcoeffhat,
            'paramhat': paramhat,
            'A_hat': A_hat,
            'B1hat': B1hat,
            'MOMvec': MOMvec,
        }

    def estimate_robustness(self, data_modifier=None, name: str = "robustness",
                            seed: int = 3991, **kwargs):
        """
        Estimate IV-VAR with data modifications for robustness checks.

        Parameters
        ----------
        data_modifier : callable that modifies X, D arrays
        name : name for this specification
        seed : random seed
        **kwargs : passed to estimate_baseline

        Returns
        -------
        Same as estimate_baseline
        """
        print(f"\n--- {name} IV-VAR Estimation ---")
        np.random.seed(seed)

        X = self.data.values[:, :self.NX].astype(np.float64)
        D = self.data.values[:, self.NX:self.NX + self.ND].astype(np.float64)

        if data_modifier is not None:
            X, D = data_modifier(X, D, **kwargs)

        MOMvec, eta_hat, A_hat = self._build_moment_vector(X, D)
        B1hat = A_hat.T.reshape(self.NX, self.NX * self.Nlags)

        param0 = self._initial_params()
        result = minimize(
            self._gmm_objective,
            param0,
            args=(MOMvec,),
            method='L-BFGS-B',
            options={'maxiter': 50000, 'ftol': 1e-15, 'gtol': 1e-10},
        )

        paramhat = result.x
        Bhat = np.zeros((self.NX, self.NX))
        Bhat[0, 0] = paramhat[0]; Bhat[1, 0] = paramhat[1]; Bhat[2, 0] = paramhat[2]
        Bhat[0, 1] = paramhat[3]; Bhat[1, 1] = paramhat[4]; Bhat[2, 1] = paramhat[5]
        Bhat[0, 2] = paramhat[6]; Bhat[1, 2] = paramhat[7]; Bhat[2, 2] = paramhat[8]

        IRF = self._compute_irf(Bhat, B1hat, X, self.lengthIRF)
        IRF_S_TO_Y = IRF[:15, 0, 2]

        print(f"  GMM objective: {result.fun:.10f}, Converged: {result.success}")

        return {
            'IRF': IRF,
            'IRF_S_TO_Y': IRF_S_TO_Y,
            'Bhat': Bhat,
            'paramhat': paramhat,
        }

    def bootstrap_se(self, baseline_result: dict, n_boot: int = 500,
                     seed: int = 3991, block_size: int = 4):
        """
        Stationary block bootstrap for standard errors.

        Matches BOOT/VAR.m with stationaryBB.m.

        Parameters
        ----------
        baseline_result : output from estimate_baseline
        n_boot : number of bootstrap replications
        seed : random seed
        block_size : expected block size for geometric distribution

        Returns
        -------
        IRF_S_TO_Y_SE : (15,) standard errors of IRF
        """
        print(f"\n--- Bootstrap SEs ({n_boot} replications) ---")
        rng = np.random.RandomState(seed)

        X = self.data.values[:, :self.NX].astype(np.float64)
        D = self.data.values[:, self.NX:self.NX + self.ND].astype(np.float64)
        T = X.shape[0]
        Nlags = self.Nlags

        # Store bootstrap IRFs
        IRF_S_TO_Y_store = np.zeros((n_boot, 15))
        boot_bad = 0

        for boot_ct in range(n_boot):
            # Stationary block bootstrap
            Xb, Db = self._stationary_block_bootstrap(X, D, rng, block_size)

            # Estimate on bootstrap sample
            try:
                MOMvec_b, eta_hat_b, A_hat_b = self._build_moment_vector(Xb, Db)
                B1hat_b = A_hat_b.T.reshape(self.NX, self.NX * Nlags)

                param0 = self._initial_params()
                result = minimize(
                    self._gmm_objective,
                    param0,
                    args=(MOMvec_b,),
                    method='L-BFGS-B',
                    options={'maxiter': 50000, 'ftol': 1e-15, 'gtol': 1e-10},
                )

                if not result.success:
                    boot_bad += 1
                    continue

                paramhat_b = result.x
                Bhat_b = np.zeros((self.NX, self.NX))
                Bhat_b[0, 0] = paramhat_b[0]; Bhat_b[1, 0] = paramhat_b[1]; Bhat_b[2, 0] = paramhat_b[2]
                Bhat_b[0, 1] = paramhat_b[3]; Bhat_b[1, 1] = paramhat_b[4]; Bhat_b[2, 1] = paramhat_b[5]
                Bhat_b[0, 2] = paramhat_b[6]; Bhat_b[1, 2] = paramhat_b[7]; Bhat_b[2, 2] = paramhat_b[8]

                # Scale factor from baseline
                SCALEFACT = (np.sqrt(np.var(X[:, 2])) /
                             np.sqrt(np.var(Xb[:, 2])) *
                             baseline_result['Bhat'][2, 2] / Bhat_b[2, 2])

                IRF_b = self._compute_irf(Bhat_b, B1hat_b, Xb, self.lengthIRF)
                IRF_b *= SCALEFACT
                IRF_S_TO_Y_store[boot_ct] = IRF_b[:15, 0, 2]

            except Exception as e:
                boot_bad += 1
                continue

            if (boot_ct + 1) % 100 == 0:
                print(f"  Bootstrap {boot_ct + 1}/{n_boot}")

        # Remove failed bootstraps
        valid = IRF_S_TO_Y_store.any(axis=1)
        IRF_S_TO_Y_store = IRF_S_TO_Y_store[valid]
        print(f"  Failed bootstraps: {boot_bad}")

        # Compute SEs
        IRF_S_TO_Y_SE = np.std(IRF_S_TO_Y_store, axis=0, ddof=1)
        print(f"  IRF_S_TO_Y_SE (first 5): {IRF_S_TO_Y_SE[:5]}")

        return IRF_S_TO_Y_SE

    def _stationary_block_bootstrap(self, X: np.ndarray, D: np.ndarray,
                                     rng: np.random.RandomState,
                                     block_size: int):
        """
        Stationary block bootstrap with geometric block sizes.

        Matches MATLAB's stationaryBB.m with sim=1 (geometric pdf).
        """
        T, NX = X.shape
        ND = D.shape[1]

        # Append first n-1 observations (wrap-around)
        X_ext = np.vstack([X, X[:T - 1]])
        D_ext = np.vstack([D, D[:T - 1]])
        n_ext = T - 1  # length of extended series

        # Generate random starting indices
        I = np.round(rng.uniform(0, 1, T) * (n_ext - 1)).astype(int)

        # Generate geometric block sizes
        b = rng.geometric(1.0 / block_size, T)

        # Bootstrap
        Xb = np.zeros((T, NX))
        Db = np.zeros((T, ND))

        for j in range(NX):
            h = 0
            for m in range(T):
                for jj in range(b[m]):
                    if h >= T:
                        break
                    Xb[h, j] = X_ext[I[m] + jj, j]
                    h += 1
                if h >= T:
                    break

        for j in range(ND):
            h = 0
            for m in range(T):
                for jj in range(b[m]):
                    if h >= T:
                        break
                    Db[h, j] = D_ext[I[m] + jj, j]
                    h += 1
                if h >= T:
                    break

        return Xb, Db

    def plot_figure6(self, IRF_BASE: np.ndarray, IRF_SE: np.ndarray):
        """
        FIGURE 6: Baseline IRF with confidence bands.

        Matches STEP2_GRAPHS.m Figure 6.
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        IRFsamp = np.arange(1, 16)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(IRFsamp, IRF_BASE, 'bo-', linewidth=2.5, label='Baseline')
        ax.plot(IRFsamp, IRF_BASE - 1.645 * IRF_SE, 'b--', linewidth=2.5)
        ax.plot(IRFsamp, IRF_BASE + 1.645 * IRF_SE, 'b--', linewidth=2.5)
        ax.plot(IRFsamp, np.zeros(15), 'k--', linewidth=2.5)
        ax.set_ylim(-10, 5)
        ax.set_xlabel('Quarters, Shock in Period 1', fontsize=14)
        ax.set_ylabel('GDP Growth, Percent Year-on-Year', fontsize=14)
        ax.tick_params(labelsize=14)

        out_path = self.output_dir / "FIGURE6.pdf"
        fig.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {out_path}")
        return fig

    def plot_figure7(self, results: dict):
        """
        FIGURE 7: Robustness IRFs.

        Matches STEP2_GRAPHS.m Figure 7.
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        IRFsamp = np.arange(1, 16)
        IRF_BASE = results['BASELINE']['IRF_S_TO_Y']
        IRF_SE = results['BOOT_SE']

        fig, ax = plt.subplots(figsize=(10, 6))

        # Confidence bands
        ax.plot(IRFsamp, IRF_BASE - 1.645 * IRF_SE, 'b--', linewidth=2.5)
        ax.plot(IRFsamp, IRF_BASE + 1.645 * IRF_SE, 'b--', linewidth=2.5)

        # Robustness specs
        specs = [
            ('FEWER_LAGS', 'c+-', 2.5),
            ('MORE_LAGS', 'gs-', 2.5),
            ('NO_TIME_FE', 'm>-', 2.5),
            ('NO_COUNTRY_FE', 'rx-', 2.5),
            ('EARLY', 'kh-', 2.5),
            ('LATE', 'yp-', 2.5),
        ]

        for name, marker, lw in specs:
            if name in results:
                ax.plot(IRFsamp, results[name]['IRF_S_TO_Y'], marker, linewidth=lw, label=name)

        # Baseline on top
        ax.plot(IRFsamp, IRF_BASE, 'bo-', linewidth=2.5, label='Baseline')
        ax.plot(IRFsamp, np.zeros(15), 'k--', linewidth=2.5)
        ax.set_ylim(-10, 5)
        ax.set_xlabel('Quarters, Shock in Period 1', fontsize=14)
        ax.set_ylabel('GDP Growth, Percent Year-on-Year', fontsize=14)
        ax.tick_params(labelsize=14)

        out_path = self.output_dir / "FIGURE7.pdf"
        fig.savefig(out_path, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved {out_path}")
        return fig

    def run_all(self):
        """Run full IV-VAR estimation pipeline."""
        self.load_data()

        # Step 1: Baseline estimation
        baseline = self.estimate_baseline()

        # Step 2: Bootstrap SEs
        boot_se = self.bootstrap_se(baseline, n_boot=500, seed=3991, block_size=4)

        # Step 3: Robustness checks
        results = {'BASELINE': baseline, 'BOOT_SE': boot_se}

        # Early/Late sample splits would require data subsetting
        # For now, run with modified lag structures
        for nlags, name in [(2, 'FEWER_LAGS'), (4, 'MORE_LAGS')]:
            orig_nlags = self.Nlags
            self.Nlags = nlags
            res = self.estimate_robustness(name=name, seed=3991)
            results[name] = res
            self.Nlags = orig_nlags

        # Step 4: Generate figures
        self.plot_figure6(baseline['IRF_S_TO_Y'], boot_se)
        self.plot_figure7(results)

        print("\n" + "=" * 70)
        print("IV-VAR ESTIMATION COMPLETE")
        print("=" * 70)

        return results
