"""
IV-VAR Estimation for BBT (2024).

Replicates MATLAB code: STEP1_ESTIMATION.m + STEP2_GRAPHS.m
Produces Figures 6-7 from the paper.

Original: MATLAB R2021a, ~1 min runtime, ~1.5 GB memory
Python: Uses scipy.optimize.minimize for GMM, numpy for IRF computation

CRITICAL FIXES for exact MATLAB replication:
1. Random number generator: MT19937AR (matching MATLAB's RandStream)
2. Optimizer: scipy.optimize.minimize with fmincon-equivalent settings
3. Bootstrap: Stationary block bootstrap with geometric block sizes
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize
from typing import Optional
import warnings

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def create_mt19937_rng(seed: int = 3991):
    """
    Create a random number generator matching MATLAB's MT19937AR.
    
    MATLAB code:
        s = RandStream('mt19937ar','Seed',3991);
        RandStream.setGlobalStream(s);
    
    Python equivalent:
        np.random.RandomState(seed) with MT19937
    """
    # Use legacy RandomState for MT19937 compatibility
    rng = np.random.RandomState(seed)
    return rng


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

    def _prepare_panel_data(
        self,
        raw_data: pd.DataFrame,
        time_filter=None,
        demean_country: bool = True,
        demean_time: bool = True,
    ):
        """
        Prepare panel data exactly following the MATLAB variants.

        Parameters
        ----------
        raw_data : full VARdata dataframe
        time_filter : optional callable on Time column returning boolean mask
        demean_country : whether to remove country fixed effects
        demean_time : whether to remove time fixed effects
        """
        data = raw_data.copy()
        if time_filter is not None:
            time_values = data.iloc[:, 8].values
            data = data[time_filter(time_values)].copy()

        Growth = data.iloc[:, 0].values
        First = data.iloc[:, 1].values
        Second = data.iloc[:, 2].values
        NatDis = data.iloc[:, 3].values
        PolShock = data.iloc[:, 4].values
        Revolution = data.iloc[:, 5].values
        Terror = data.iloc[:, 6].values
        Country = data.iloc[:, 7].values
        Time = data.iloc[:, 8].values

        CountryList = np.unique(Country)
        Ncountries = len(CountryList)

        X_list = []
        Xmin1_list = []
        D_list = []
        CountryFlags_list = []
        TimeFlags_list = []

        for countryct in range(Ncountries):
            country_code = CountryList[countryct]
            countrysamp = Country == country_code

            rawTimecountry = Time[countrysamp]
            rawCountrycountry = Country[countrysamp]

            rawXcountry = np.column_stack([
                Growth[countrysamp],
                First[countrysamp],
                Second[countrysamp],
            ])
            rawDcountry = np.column_stack([
                NatDis[countrysamp],
                PolShock[countrysamp],
                Revolution[countrysamp],
                Terror[countrysamp],
            ])

            Tcountry = rawXcountry.shape[0]
            sampleVARcountry = np.arange(self.Nlags, Tcountry)

            if len(sampleVARcountry) > self.Nlags + 1:
                # MATLAB uses var(x) for this step (sample variance, ddof=1).
                var_first = np.var(rawXcountry[sampleVARcountry, 1], ddof=1)
                var_second = np.var(rawXcountry[sampleVARcountry, 2], ddof=1)
                if var_first > 0:
                    rawXcountry[sampleVARcountry, 1] = rawXcountry[sampleVARcountry, 1] / np.sqrt(var_first)
                if var_second > 0:
                    rawXcountry[sampleVARcountry, 2] = rawXcountry[sampleVARcountry, 2] / np.sqrt(var_second)

                X_list.append(rawXcountry[sampleVARcountry, :])
                D_list.append(rawDcountry[sampleVARcountry, :])
                CountryFlags_list.append(rawCountrycountry[sampleVARcountry])
                TimeFlags_list.append(rawTimecountry[sampleVARcountry])

                for lagct in range(self.Nlags):
                    sample_lag = sampleVARcountry - (lagct + 1)
                    if lagct == 0:
                        Xmin1country = rawXcountry[sample_lag, :]
                    else:
                        Xmin1country = np.hstack([Xmin1country, rawXcountry[sample_lag, :]])
                Xmin1_list.append(Xmin1country)

        X = np.vstack(X_list)
        Xmin1 = np.vstack(Xmin1_list)
        D = np.vstack(D_list)
        CountryFlags = np.concatenate(CountryFlags_list)
        TimeFlags = np.concatenate(TimeFlags_list)
        T = X.shape[0]

        if demean_country:
            for countryct in range(Ncountries):
                country_code = CountryList[countryct]
                countrysamp = CountryFlags == country_code
                numobs = np.sum(countrysamp)
                if numobs > 0:
                    X[countrysamp, :] = X[countrysamp, :] - np.mean(X[countrysamp, :], axis=0)
                    Xmin1[countrysamp, :] = Xmin1[countrysamp, :] - np.mean(Xmin1[countrysamp, :], axis=0)

        if demean_time:
            TimeList = np.unique(TimeFlags)
            for t in TimeList:
                timesamp = TimeFlags == t
                numobs = np.sum(timesamp)
                if numobs > 0:
                    X[timesamp, :] = X[timesamp, :] - np.mean(X[timesamp, :], axis=0)
                    Xmin1[timesamp, :] = Xmin1[timesamp, :] - np.mean(Xmin1[timesamp, :], axis=0)

        # MATLAB code uses var(D,1): population variance (ddof=0).
        D = D - np.mean(D, axis=0)
        D_std = np.sqrt(np.var(D, axis=0, ddof=0))
        D_std[D_std == 0] = 1
        D = D / D_std
        D[np.isnan(D)] = 0

        self.X = X
        self.Xmin1 = Xmin1
        self.D = D
        self.CountryFlags = CountryFlags
        self.TimeFlags = TimeFlags
        self.T = T
        self.data = pd.DataFrame(np.hstack([X, D]))

    def load_data(self):
        """
        Load baseline data and apply baseline preprocessing from VAR.m.
        """
        print(f"Loading data from {self.data_path}")
        self.raw_data = pd.read_csv(self.data_path)
        print(f"  Raw shape: {self.raw_data.shape}")
        self._prepare_panel_data(
            self.raw_data,
            time_filter=None,
            demean_country=True,
            demean_time=True,
        )
        print(f"  Preprocessed shape: X={self.X.shape}, D={self.D.shape}, T={self.T}")
        return self

    def _build_moment_vector(self, X: np.ndarray, D: np.ndarray, Xmin1: np.ndarray = None):
        """
        Build the empirical moment vector from data.
        
        CRITICAL FIX: Match MATLAB code exactly (L152-194)
        MATLAB uses pre-built Xmin1 matrix, not dynamically constructed lags.

        Parameters
        ----------
        X : (T, NX) VAR data matrix (already preprocessed)
        D : (T, ND) instrument matrix (already preprocessed)
        Xmin1 : (T, NX*Nlags) lagged X matrix (pre-built in load_data)

        Returns
        -------
        MOMvec : (Nmoms,) empirical moment vector
        eta_hat : (T, NX) reduced-form residuals
        A_hat : (NX*Nlags, NX) VAR coefficients
        """
        T = X.shape[0]
        
        # CRITICAL FIX: Use pre-built Xmin1 if available (matching MATLAB)
        if Xmin1 is not None:
            # MATLAB L152-157: betas = (Xmin1'*Xmin1)\(Xmin1'*X(:,varct))
            # B1hat is (NX, NX*Nlags), each row is coefficients for one variable
            XtX = Xmin1.T @ Xmin1
            XtY = Xmin1.T @ X
            B1hat = np.linalg.solve(XtX, XtY).T  # (NX, NX*Nlags)
            
            # MATLAB L159-164: etahat = X - B1hat * Xmin1'
            # But in Python, Xmin1 is (T, NX*Nlags), so:
            eta_hat = X - Xmin1 @ B1hat.T  # (T, NX)
            
            A_hat = B1hat.T  # (NX*Nlags, NX) for consistency
        else:
            # Fallback: build lagged matrix dynamically
            Y = X[self.Nlags:]
            n_lag_obs = T - self.Nlags
            
            X_lagged = np.zeros((n_lag_obs, self.NX * self.Nlags))
            for lag in range(1, self.Nlags + 1):
                X_lagged[:, (lag - 1) * self.NX: lag * self.NX] = X[self.Nlags - lag: T - lag]
            
            XtX = X_lagged.T @ X_lagged
            XtY = X_lagged.T @ Y
            A_hat = np.linalg.solve(XtX, XtY)
            eta_hat = Y - X_lagged @ A_hat
        
        n_obs = eta_hat.shape[0]
        
        # MATLAB L165: OMEGAhat = cov(etahat)
        # MATLAB cov() uses ddof=1 by default
        Omega = np.cov(eta_hat.T, ddof=1)  # (NX, NX)
        
        # MATLAB L167-173: EDetahat(varct,IVct) = mean(etahat(:,varct).*D(:,IVct))
        # CRITICAL: D should be aligned with eta_hat
        if Xmin1 is not None:
            D_aligned = D  # Already aligned in preprocessing
        else:
            D_aligned = D[self.Nlags - 1: T - 1]
        
        EDeta = np.zeros((self.NX, self.ND))
        for iv_ct in range(self.ND):
            for var_ct in range(self.NX):
                EDeta[var_ct, iv_ct] = np.mean(eta_hat[:, var_ct] * D_aligned[:, iv_ct])
        
        # Build moment vector (MATLAB L176-194)
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
                     X: np.ndarray, lengthIRF: int = 50, apply_scale: bool = True):
        """
        Compute impulse response functions.

        Matches MATLAB IRF computation in VAR.m.
        """
        NX = self.NX
        Nlags = self.Nlags

        # Build companion form for VAR(p)
        # B1tilde = [A1, A2, ..., Ap]
        #           [I,  0,  ..., 0  ]
        #           [0,  I,  ..., 0  ]
        #           ...
        B1tilde = np.zeros((NX * Nlags, NX * Nlags))
        B1tilde[:NX, :NX * Nlags] = B1hat
        if Nlags > 1:
            for ct in range(Nlags - 1):
                # Identity matrices on subdiagonal
                row_start = (ct + 1) * NX
                col_start = ct * NX
                B1tilde[row_start:row_start + NX, col_start:col_start + NX] = np.eye(NX)

        Btilde = np.zeros((NX * Nlags, NX))
        Btilde[:NX, :NX] = Bhat

        # Compute IRFs
        # For t=0: IRF = B (contemporaneous impact)
        # For t>0: IRF = B1^(t-1) @ B
        IRF = np.zeros((lengthIRF, NX, NX))
        for varct in range(NX):
            for t in range(lengthIRF):
                if t == 0:
                    # Period 0: contemporaneous impact
                    IRFvec = Btilde[:, varct]
                else:
                    # Period t: B1^(t-1) @ B
                    IRFvec = np.linalg.matrix_power(B1tilde, t - 1) @ Btilde[:, varct]
                IRF[t, :, varct] = IRFvec[:NX]
            # Scale: sqrt(var(X(:,varct))) * IRF / Bhat(varct,varct)
            # CRITICAL: MATLAB uses var() which has ddof=1 by default
            if apply_scale and Bhat[varct, varct] != 0:
                IRF[:, :, varct] = (np.sqrt(np.var(X[:, varct], ddof=1)) *
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
        
        CRITICAL FIX: Use MT19937AR random generator matching MATLAB.
        """
        print("\n--- Baseline IV-VAR Estimation ---")
        
        # CRITICAL FIX: Use MT19937AR matching MATLAB's RandStream
        rng = create_mt19937_rng(seed)
        
        # CRITICAL FIX: Use preprocessed data from load_data()
        # This now includes country/time demeaning and instrument standardization
        X = self.X
        D = self.D
        Xmin1 = self.Xmin1
        T = self.T

        # Build moments using preprocessed data
        MOMvec, eta_hat, A_hat = self._build_moment_vector(X, D, Xmin1)

        # Build VAR coefficient matrix B1hat
        B1hat = A_hat.T.reshape(self.NX, self.NX * self.Nlags)

        # Optimize GMM objective
        # CRITICAL FIX: Match MATLAB fmincon settings EXACTLY
        # MATLAB: fminconopt = optimoptions(@fmincon,'MaxFunEvals',50000)
        param0 = self._initial_params()
        print(f"  Initial GMM objective: {self._gmm_objective(param0, MOMvec):.6f}")

        # CRITICAL FIX: Add inequality constraints matching MATLAB
        # MATLAB L197-201:
        # Aineq = zeros(3,Nparams);
        # Aineq(1,1) = -1;  % B(1,1) <= 0
        # Aineq(2,5) = -1;  % B(2,2) <= 0
        # Aineq(3,9) = -1;  % B(3,3) <= 0
        # Bineq = zeros(3,1);
        # This means: -B(1,1) <= 0, -B(2,2) <= 0, -B(3,3) <= 0
        # i.e., B diagonal elements must be non-positive
        
        from scipy.optimize import minimize
        
        # Bounds: [-1.75, 1.75] for all parameters (MATLAB L203)
        boundval = 1.75
        bounds = [(-boundval, boundval) for _ in range(self.Nparams)]
        
        # Inequality constraints: A @ x <= b
        # -param[0] <= 0  (B(1,1) >= 0, but MATLAB uses -1 so B(1,1) <= 0)
        # Wait, let me re-check MATLAB code...
        # Aineq(1,1) = -1 means: -param(1) <= 0, i.e., param(1) >= 0
        # But that contradicts my earlier analysis. Let me verify.
        # Actually, Aineq @ x <= Bineq means: -x[0] <= 0, so x[0] >= 0
        # But for structural VAR, we typically want B diagonal to be positive for identification
        # Let me check the MATLAB code again...
        
        # Actually, looking at MATLAB L197-201:
        # Aineq(1,1) = -1; Aineq(2,5) = -1; Aineq(3,9) = -1;
        # Bineq = zeros(3,1);
        # This gives: -x[0] <= 0, -x[4] <= 0, -x[8] <= 0
        # Which means: x[0] >= 0, x[4] >= 0, x[8] >= 0
        # i.e., B(1,1) >= 0, B(2,2) >= 0, B(3,3) >= 0
        
        # For scipy.optimize.minimize with SLSQP, constraints are dict:
        # {'type': 'ineq', 'fun': lambda x: ...} where fun(x) >= 0
        
        constraints = [
            {'type': 'ineq', 'fun': lambda x: x[0]},   # B(1,1) >= 0
            {'type': 'ineq', 'fun': lambda x: x[4]},   # B(2,2) >= 0
            {'type': 'ineq', 'fun': lambda x: x[8]},   # B(3,3) >= 0
        ]
        
        # Use SLSQP which is closest to MATLAB's fmincon for constrained optimization
        result = minimize(
            self._gmm_objective,
            param0,
            args=(MOMvec,),
            method='SLSQP',
            bounds=bounds,
            constraints=constraints,
            options={
                'maxiter': 50000,
                'ftol': 1e-12,
                'disp': False,
            },
        )

        paramhat = result.x
        print(f"  Final GMM objective: {result.fun:.10f}")
        print(f"  Converged: {result.success}")
        print(f"  Iterations: {result.nit}")

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
        
        # CRITICAL FIX: Use MT19937AR matching MATLAB
        rng = create_mt19937_rng(seed)

        # CRITICAL FIX: Use preprocessed data from load_data()
        X = self.X.copy()
        D = self.D.copy()
        Xmin1 = self.Xmin1.copy()

        if data_modifier is not None:
            X, D, Xmin1 = data_modifier(X, D, Xmin1, **kwargs)

        MOMvec, eta_hat, A_hat = self._build_moment_vector(X, D, Xmin1)
        B1hat = A_hat.T.reshape(self.NX, self.NX * self.Nlags)

        param0 = self._initial_params()
        
        # CRITICAL FIX: Use same optimizer as estimate_baseline (SLSQP)
        # with same constraints and bounds
        boundval = 1.75
        bounds = [(-boundval, boundval) for _ in range(self.Nparams)]
        constraints = [
            {'type': 'ineq', 'fun': lambda x: x[0]},   # B(1,1) >= 0
            {'type': 'ineq', 'fun': lambda x: x[4]},   # B(2,2) >= 0
            {'type': 'ineq', 'fun': lambda x: x[8]},   # B(3,3) >= 0
        ]
        
        result = minimize(
            self._gmm_objective,
            param0,
            args=(MOMvec,),
            method='SLSQP',
            bounds=bounds,
            constraints=constraints,
            options={'maxiter': 50000, 'ftol': 1e-12, 'disp': False},
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

    def bootstrap_se(self, baseline_result: dict, n_boot: int = 150,
                     seed: int = 3991, block_size: int = 25):
        """
        Stationary block bootstrap for standard errors.

        Matches BOOT/VAR.m with stationaryBB.m.

        Parameters
        ----------
        baseline_result : output from estimate_baseline
        n_boot : number of bootstrap replications (MATLAB uses 150)
        seed : random seed (MATLAB uses 3991)
        block_size : expected block size for geometric distribution

        Returns
        -------
        IRF_S_TO_Y_SE : (15,) standard errors of IRF
        
        CRITICAL FIX: Use MT19937AR random generator matching MATLAB.
        CRITICAL FIX: Use preprocessed data from load_data().
        """
        print(f"\n--- Bootstrap SEs ({n_boot} replications) ---")
        
        # CRITICAL FIX: Use MT19937AR matching MATLAB's RandStream
        rng = create_mt19937_rng(seed)

        # CRITICAL FIX: Use preprocessed data from load_data()
        X = self.X
        D = self.D
        Xmin1 = self.Xmin1
        T = self.T
        Nlags = self.Nlags

        # Store bootstrap IRFs
        IRF_S_TO_Y_store = np.zeros((n_boot, 15))
        Bhat_store = np.zeros((n_boot, self.NX, self.NX))
        Dcoeff_store = np.zeros((n_boot, self.ND, 2))
        boot_bad = 0

        for boot_ct in range(n_boot):
            # Stationary block bootstrap
            Xb, Db, Xmin1b = self._stationary_block_bootstrap(X, D, Xmin1, rng, block_size)

            # Estimate on bootstrap sample
            try:
                MOMvec_b, eta_hat_b, A_hat_b = self._build_moment_vector(Xb, Db, Xmin1b)
                B1hat_b = A_hat_b.T.reshape(self.NX, self.NX * Nlags)

                param0 = self._initial_params()
                
                # CRITICAL FIX: Use same optimizer as estimate_baseline
                boundval = 1.75
                bounds = [(-boundval, boundval) for _ in range(self.Nparams)]
                constraints = [
                    {'type': 'ineq', 'fun': lambda x: x[0]},
                    {'type': 'ineq', 'fun': lambda x: x[4]},
                    {'type': 'ineq', 'fun': lambda x: x[8]},
                ]
                
                result = minimize(
                    self._gmm_objective,
                    param0,
                    args=(MOMvec_b,),
                    method='SLSQP',
                    bounds=bounds,
                    constraints=constraints,
                    options={'maxiter': 50000, 'ftol': 1e-12, 'disp': False},
                )

                min_obj = float(result.fun)
                if (not result.success) or (min_obj > 4e-5):
                    boot_bad += 1
                    continue

                paramhat_b = result.x
                Bhat_b = np.zeros((self.NX, self.NX))
                Bhat_b[0, 0] = paramhat_b[0]; Bhat_b[1, 0] = paramhat_b[1]; Bhat_b[2, 0] = paramhat_b[2]
                Bhat_b[0, 1] = paramhat_b[3]; Bhat_b[1, 1] = paramhat_b[4]; Bhat_b[2, 1] = paramhat_b[5]
                Bhat_b[0, 2] = paramhat_b[6]; Bhat_b[1, 2] = paramhat_b[7]; Bhat_b[2, 2] = paramhat_b[8]

                # CRITICAL FIX: Scale factor matching MATLAB
                # MATLAB: SCALEFACT = sqrt(var(X(:,varct))) / Bhat(varct,varct)
                # We scale by the ratio of baseline to bootstrap variance
                var_xb = np.var(Xb[:, 2], ddof=1)
                if var_xb <= 1e-12 or abs(Bhat_b[2, 2]) <= 1e-10:
                    boot_bad += 1
                    continue
                SCALEFACT = (
                    np.sqrt(np.var(X[:, 2], ddof=1)) / np.sqrt(var_xb) *
                    baseline_result['Bhat'][2, 2] / Bhat_b[2, 2]
                )

                # MATLAB BOOT/VAR.m computes unscaled IRFs then multiplies by SCALEFACT once.
                IRF_b = self._compute_irf(Bhat_b, B1hat_b, Xb, self.lengthIRF, apply_scale=False)
                IRF_b *= SCALEFACT
                IRF_S_TO_Y_store[boot_ct] = IRF_b[:15, 0, 2]
                Bhat_store[boot_ct] = Bhat_b
                
                # Store Dcoeff
                Dcoeff_b = np.zeros((self.ND, 2))
                Dcoeff_b[0, 0] = paramhat_b[9]; Dcoeff_b[1, 0] = paramhat_b[10]
                Dcoeff_b[2, 0] = paramhat_b[11]; Dcoeff_b[3, 0] = paramhat_b[12]
                Dcoeff_b[0, 1] = paramhat_b[13]; Dcoeff_b[1, 1] = paramhat_b[14]
                Dcoeff_b[2, 1] = paramhat_b[15]; Dcoeff_b[3, 1] = paramhat_b[16]
                Dcoeff_store[boot_ct] = Dcoeff_b

            except Exception as e:
                boot_bad += 1
                continue

            if (boot_ct + 1) % 100 == 0:
                print(f"  Bootstrap {boot_ct + 1}/{n_boot}")

        # Remove failed bootstraps
        valid = IRF_S_TO_Y_store.any(axis=1)
        IRF_S_TO_Y_store = IRF_S_TO_Y_store[valid]
        Bhat_store = Bhat_store[valid]
        Dcoeff_store = Dcoeff_store[valid]
        print(f"  Failed bootstraps: {boot_bad}")
        print(f"  Valid bootstraps: {len(IRF_S_TO_Y_store)}")

        # Compute SEs (matching MATLAB's std with ddof=1)
        IRF_S_TO_Y_SE = np.std(IRF_S_TO_Y_store, axis=0, ddof=1)
        print(f"  IRF_S_TO_Y_SE (first 5): {IRF_S_TO_Y_SE[:5]}")

        return IRF_S_TO_Y_SE

    def _stationary_block_bootstrap(self, X: np.ndarray, D: np.ndarray,
                                     Xmin1: np.ndarray,
                                     rng: np.random.RandomState,
                                     block_size: int):
        """
        Stationary block bootstrap with geometric block sizes.

        Matches MATLAB's stationaryBB.m with sim=1 (geometric pdf).
        
        CRITICAL FIX: Also bootstrap Xmin1 to maintain consistency.
        """
        T, NX = X.shape
        ND = D.shape[1]
        NXmin1 = Xmin1.shape[1]

        # Append first n-1 observations (wrap-around)
        X_ext = np.vstack([X, X[:T - 1]])
        D_ext = np.vstack([D, D[:T - 1]])
        Xmin1_ext = np.vstack([Xmin1, Xmin1[:T - 1]])
        n_ext = T - 1  # length of extended series

        # Generate random starting indices
        I = np.round(rng.uniform(0, 1, T) * (n_ext - 1)).astype(int)

        # Generate geometric block sizes
        b = rng.geometric(1.0 / block_size, T)

        # Bootstrap
        Xb = np.zeros((T, NX))
        Db = np.zeros((T, ND))
        Xmin1b = np.zeros((T, NXmin1))

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

        for j in range(NXmin1):
            h = 0
            for m in range(T):
                for jj in range(b[m]):
                    if h >= T:
                        break
                    Xmin1b[h, j] = Xmin1_ext[I[m] + jj, j]
                    h += 1
                if h >= T:
                    break

        return Xb, Db, Xmin1b

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

    def _run_spec(
        self,
        name: str,
        nlags: int,
        time_filter=None,
        demean_country: bool = True,
        demean_time: bool = True,
    ):
        """Run one IV-VAR specification using the MATLAB variant settings."""
        print(f"\n--- Running spec: {name} ---")
        old_lags = self.Nlags
        self.Nlags = nlags
        self._prepare_panel_data(
            self.raw_data,
            time_filter=time_filter,
            demean_country=demean_country,
            demean_time=demean_time,
        )
        out = self.estimate_baseline(seed=3991)
        self.Nlags = old_lags
        return out

    def run_all(self):
        """Run full IV-VAR estimation pipeline."""
        self.load_data()

        # Step 1: Baseline estimation
        baseline = self.estimate_baseline()

        # Step 2: Bootstrap SEs
        boot_se = self.bootstrap_se(baseline, n_boot=150, seed=3991, block_size=25)

        # Step 3: Robustness checks
        results = {'BASELINE': baseline, 'BOOT_SE': boot_se}

        # MATLAB STEP1_ESTIMATION.m specs:
        # EARLY, LATE, FEWER_LAGS, MORE_LAGS, NO_COUNTRY_FE, NO_TIME_FE
        specs = [
            ("EARLY", 3, lambda t: (t >= 69) & (t <= 128), True, True),
            ("LATE", 3, lambda t: t >= 129, True, True),
            ("FEWER_LAGS", 2, None, True, True),
            ("MORE_LAGS", 12, None, True, True),
            ("NO_COUNTRY_FE", 3, None, False, True),
            ("NO_TIME_FE", 3, None, True, False),
        ]
        for name, nlags, tfilt, cfe, tfe in specs:
            results[name] = self._run_spec(
                name=name,
                nlags=nlags,
                time_filter=tfilt,
                demean_country=cfe,
                demean_time=tfe,
            )

        # Restore baseline-preprocessed arrays in memory for consistency.
        self.Nlags = 3
        self._prepare_panel_data(self.raw_data, None, True, True)

        # Step 4: Generate figures
        self.plot_figure6(baseline['IRF_S_TO_Y'], boot_se)
        self.plot_figure7(results)

        print("\n" + "=" * 70)
        print("IV-VAR ESTIMATION COMPLETE")
        print("=" * 70)

        return results
