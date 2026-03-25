"""
Micro-Macro Model Simulation for BBT (2024).

Replicates Fortran 90 code: VOL_GROWTH_wrapper.f90 + base_lib.f90
Produces Figure 8 from the paper.

The model is a heterogeneous-firm DSGE with:
- Firms facing time-varying uncertainty (volatility shocks)
- NON-CONVEX adjustment costs (fixed costs + linear costs)
- Aggregate dynamics emerging from micro-level decisions

Key Features:
1. NON-CONVEX adjustment costs for capital and labor
2. 5D state space: (z, a, s, k, l)
3. Howard acceleration for VFI
4. GMM estimation with PSO optimization
5. IV regression for moment matching

Original: Fortran 90, ~5 min runtime
Python: Uses numpy/scipy for grid solving
"""

import sys
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional
import warnings

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from .params import ModelParameters, create_params
from .grids import StateGrids, build_grids
from .adjustment import (
    AdjustmentCostCalculator,
    capital_adjustment_cost,
    labor_adjustment_cost,
    output
)
from .vfi import VFISolution, solve_vfi, solve_vfi_simplified
from .simulation import (
    SimulationResults,
    simulate_firms,
    simulate_irf,
    compute_figure8_irf
)
from .gmm import GMMSolution, estimate_gmm, gmm_objective
from .iv_regression import IVResults, run_iv_regression
from .optimizer import PSOConfig, PSOResult, pso_optimize


class MicroMacroModel:
    """
    Micro-macro model matching the Fortran implementation.
    
    This class provides a unified interface to:
    1. Build state space grids
    2. Solve firm value function (VFI)
    3. Simulate firm dynamics
    4. Compute IRFs
    5. Estimate parameters via GMM
    
    Examples
    --------
    >>> model = MicroMacroModel(simplified=True)
    >>> model.run_all()
    """
    
    def __init__(self, simplified: bool = True):
        """
        Initialize the model.
        
        Parameters
        ----------
        simplified : bool
            If True, use smaller grid sizes for faster computation.
        """
        self.output_dir = PROJECT_ROOT / "output" / "figures"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.simplified = simplified
        
        # Initialize parameters
        self.params = create_params(simplified=simplified)
        
        # Initialize grids
        self.grids: Optional[StateGrids] = None
        
        # VFI solution
        self.vfi_solution: Optional[VFISolution] = None
        
        # Simulation results
        self.sim_results: Optional[SimulationResults] = None
        self.irf_results = None
        
        # GMM results
        self.gmm_solution: Optional[GMMSolution] = None
        
    def build(self):
        """Build state space grids and pre-compute matrices."""
        print("\n--- Building State Space ---")
        self.grids = build_grids(self.params)
        print(f"  Grid sizes: z={self.params.znum}, a={self.params.anum}, "
              f"s={self.params.snum}, k={self.params.knum}, l={self.params.lnum}")
        print(f"  Total states: {self.grids.numstates:,}")
        
    def solve(self, max_iter: int = None, tol: float = None, verbose: bool = True):
        """
        Solve firm value function via VFI.
        
        Parameters
        ----------
        max_iter : int, optional
            Max VFI iterations.
        tol : float, optional
            Convergence tolerance.
        verbose : bool
            Print progress.
        """
        if self.grids is None:
            self.build()
        
        print("\n--- Solving Value Function ---")
        
        # Use simplified VFI for speed
        self.vfi_solution = solve_vfi_simplified(
            self.params, self.grids,
            max_iter=max_iter or self.params.vfmaxit,
            tol=tol or self.params.vferrortol,
            verbose=verbose
        )
        
        print(f"  Converged: {self.vfi_solution.converged}")
        print(f"  Iterations: {self.vfi_solution.iterations}")
        print(f"  VF error: {self.vfi_solution.vf_error:.2e}")
        
    def simulate(self, T: int = None, seed: int = 2501):
        """
        Simulate firm dynamics and aggregate variables.
        
        Parameters
        ----------
        T : int, optional
            Number of periods.
        seed : int
            Random seed.
        """
        if self.vfi_solution is None:
            self.solve()
        
        print("\n--- Running Simulation ---")
        
        self.sim_results = simulate_firms(
            self.params, self.grids, self.vfi_solution,
            T=T, seed=seed
        )
        
        print(f"  Simulated {len(self.sim_results.Y_sim)} periods")
        print(f"  Mean GDP: {np.mean(self.sim_results.Y_sim):.4f}")
        print(f"  Mean Investment: {np.mean(self.sim_results.I_sim):.4f}")
        
    def compute_irf(self, T: int = 40, n_sims: int = 100):
        """
        Compute impulse response to uncertainty shock.
        
        Parameters
        ----------
        T : int
            IRF length in periods.
        n_sims : int
            Number of simulations to average.
        """
        if self.vfi_solution is None:
            self.solve()
        
        print("\n--- Computing IRF ---")
        
        self.irf_results = simulate_irf(
            self.params, self.grids, self.vfi_solution,
            T_irf=T, n_sims=n_sims
        )
        
        print(f"  Max GDP drop: {np.min(self.irf_results.irf_Y):.2f}%")
        print(f"  Max Investment drop: {np.min(self.irf_results.irf_I):.2f}%")
        
    def estimate(self, max_evals: int = 100):
        """
        Estimate model parameters via GMM.
        
        Parameters
        ----------
        max_evals : int
            Maximum function evaluations.
        """
        if self.grids is None:
            self.build()
        
        print("\n--- GMM Estimation ---")
        
        self.gmm_solution = estimate_gmm(
            self.params, self.grids,
            max_evals=max_evals,
            verbose=True
        )
        
        print(f"  Converged: {self.gmm_solution.converged}")
        print(f"  GMM value: {self.gmm_solution.gmm_value:.4f}")
        print(f"  Optimal parameters: {self.gmm_solution.x_opt}")
        
    def plot_figure8(self, save: bool = True):
        """
        Plot Figure 8: GDP and Investment IRFs.
        
        Parameters
        ----------
        save : bool
            Whether to save the figure.
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        if self.irf_results is None:
            self.compute_irf()
        
        periods = self.irf_results.periods
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # GDP response
        ax1.plot(periods, self.irf_results.irf_Y, 'bo-', linewidth=2.5)
        ax1.plot(periods, np.zeros(len(periods)), 'k--', linewidth=2.5)
        ax1.set_xlabel('Quarters', fontsize=14)
        ax1.set_ylabel('GDP, Percent Deviation', fontsize=14)
        ax1.set_title('GDP Response to Uncertainty Shock', fontsize=14)
        ax1.tick_params(labelsize=12)
        ax1.grid(True, alpha=0.3)
        
        # Investment response
        ax2.plot(periods, self.irf_results.irf_I, 'bo-', linewidth=2.5)
        ax2.plot(periods, np.zeros(len(periods)), 'k--', linewidth=2.5)
        ax2.set_xlabel('Quarters', fontsize=14)
        ax2.set_ylabel('Investment, Percent Deviation', fontsize=14)
        ax2.set_title('Investment Response to Uncertainty Shock', fontsize=14)
        ax2.tick_params(labelsize=12)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save:
            out_path = self.output_dir / "FIGURE8.pdf"
            fig.savefig(out_path, bbox_inches='tight')
            print(f"\n  Saved {out_path}")
        
        plt.close(fig)
        return fig
    
    def run_all(self, do_estimation: bool = False):
        """
        Run full model simulation pipeline.
        
        Parameters
        ----------
        do_estimation : bool
            Whether to run GMM estimation (time-consuming).
        """
        print("\n" + "=" * 70)
        print("Micro-Macro Model Simulation")
        print("=" * 70)
        
        # Print parameters
        print("\n--- Model Parameters (matching Fortran) ---")
        p = self.params
        print(f"  alpha (capital share): {p.alpha}")
        print(f"  nu (labor share): {p.nu}")
        print(f"  capirrev (capital irreversibility): {p.capirrev}")
        print(f"  hirelin (hiring cost): {p.hirelin}")
        print(f"  firelin (firing cost): {p.firelin}")
        print(f"  labfix (fixed labor adj cost): {p.labfix}")
        print(f"  beta (discount factor): {p.beta:.4f}")
        
        print("\n--- NON-CONVEX Adjustment Costs (Key Feature) ---")
        print("  Capital: ACk = capirrev * |k' - k| + capfix * I(k' != k)")
        print("  Labor: ACl = hirelin * max(l' - l, 0) + firelin * max(l - l', 0)")
        print("              + labfix * I(l' != l)")
        print("")
        print("  This is the CRITICAL difference from quadratic adjustment costs!")
        
        # Run steps
        self.build()
        self.solve()
        self.simulate()
        self.compute_irf()
        self.plot_figure8()
        
        if do_estimation:
            self.estimate()
        
        # Summary
        print("\n--- IRF Results ---")
        print(f"  Max GDP drop: {np.min(self.irf_results.irf_Y):.2f}%")
        print(f"  Max Investment drop: {np.min(self.irf_results.irf_I):.2f}%")
        
        print("\n" + "=" * 70)
        print("MODEL SIMULATION COMPLETE")
        print("=" * 70)
        
        return {
            'irf_gdp': self.irf_results.irf_Y,
            'irf_investment': self.irf_results.irf_I,
            'vfi_solution': self.vfi_solution,
            'sim_results': self.sim_results
        }


# Convenience function for quick testing
def quick_test():
    """Run a quick test of the model."""
    model = MicroMacroModel(simplified=True)
    model.solve()
    model.compute_irf(T=20, n_sims=10)
    return model


if __name__ == '__main__':
    model = MicroMacroModel(simplified=True)
    model.run_all()
