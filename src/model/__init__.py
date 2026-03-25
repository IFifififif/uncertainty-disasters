"""
MODEL Module for BBT (2024) - Using Disasters to Estimate the Impact of Uncertainty.

This module implements the heterogeneous-firm DSGE model with:
- Non-convex adjustment costs for capital and labor
- Time-varying uncertainty (volatility shocks)
- Aggregate dynamics from micro-level decisions
- GMM estimation with PSO optimization

Main Classes:
- MicroMacroModel: Main interface for running the model
- ModelParameters: All calibrated parameters
- StateGrids: State space grids and transition matrices
- VFISolution: Value function iteration results
- GMMSolution: GMM estimation results

Main Functions:
- build_grids(): Build state space
- solve_vfi(): Solve value function
- simulate_firms(): Simulate firm dynamics
- compute_figure8_irf(): Compute IRF for Figure 8
- estimate_gmm(): Estimate parameters via GMM
"""

from .params import (
    ModelParameters,
    create_params
)

from .grids import (
    StateGrids,
    build_grids,
    get_state_indices
)

from .adjustment import (
    AdjustmentCostCalculator,
    capital_adjustment_cost,
    labor_adjustment_cost,
    output,
    compute_adjustment_costs_grid,
    compute_investment_matrix,
    compute_output_matrix
)

from .vfi import (
    VFISolution,
    solve_vfi,
    solve_vfi_simplified
)

from .simulation import (
    SimulationResults,
    IRFResults,
    simulate_firms,
    simulate_irf,
    simulate_firms_with_shock,
    compute_figure8_irf
)

from .gmm import (
    GMMSolution,
    compute_simulated_moments,
    gmm_objective,
    estimate_gmm
)

from .iv_regression import (
    FirstStageResults,
    SecondStageResults,
    IVResults,
    first_stage_regression,
    second_stage_regression,
    run_iv_regression,
    compute_iv_moments,
    simulate_iv_data
)

from .optimizer import (
    PSOConfig,
    PSOResult,
    pso_optimize,
    pso_optimize_restart,
    nelder_mead_2d
)

from .solve import (
    MicroMacroModel,
    quick_test
)


__all__ = [
    # Main classes
    'MicroMacroModel',
    'ModelParameters',
    'StateGrids',
    'VFISolution',
    'GMMSolution',
    'SimulationResults',
    'IRFResults',
    'IVResults',
    'PSOConfig',
    'PSOResult',
    
    # Parameters
    'create_params',
    
    # Grids
    'build_grids',
    'get_state_indices',
    
    # Adjustment costs
    'AdjustmentCostCalculator',
    'capital_adjustment_cost',
    'labor_adjustment_cost',
    'output',
    'compute_adjustment_costs_grid',
    'compute_investment_matrix',
    'compute_output_matrix',
    
    # VFI
    'solve_vfi',
    'solve_vfi_simplified',
    
    # Simulation
    'simulate_firms',
    'simulate_irf',
    'simulate_firms_with_shock',
    'compute_figure8_irf',
    
    # GMM
    'compute_simulated_moments',
    'gmm_objective',
    'estimate_gmm',
    
    # IV Regression
    'FirstStageResults',
    'SecondStageResults',
    'first_stage_regression',
    'second_stage_regression',
    'run_iv_regression',
    'compute_iv_moments',
    'simulate_iv_data',
    
    # Optimization
    'pso_optimize',
    'pso_optimize_restart',
    'nelder_mead_2d',
    
    # Convenience
    'quick_test'
]
