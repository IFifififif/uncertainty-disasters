"""
Model Parameters for BBT (2024) - Using Disasters to Estimate the Impact of Uncertainty.

This module contains all parameters from the original Fortran VOL_GROWTH_wrapper.f90.
Parameters are organized into categories matching the original code structure.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import numpy as np


@dataclass
class ModelParameters:
    """
    Complete parameter set matching Fortran VOL_GROWTH_wrapper.f90.
    
    Parameter organization:
    1. Production function
    2. Depreciation rates
    3. Discount factor
    4. Adjustment costs (NON-CONVEX)
    5. Uncertainty process
    6. Productivity processes
    7. Grid sizes
    8. Simulation controls
    9. Disaster parameters (for GMM estimation)
    """
    
    # =====================
    # Production Function
    # =====================
    alpha: float = 0.25      # Capital share (Fortran: alpha = 0.25)
    nu: float = 0.5          # Labor share (Fortran: nu = 0.5)
    theta: float = 2.0       # Labor supply elasticity
    
    # =====================
    # Depreciation Rates
    # =====================
    deltak: float = 0.026    # Capital depreciation
    deltan: float = 0.088    # Labor depreciation (worker separations)
    
    # =====================
    # Discount Factor
    # =====================
    # Fortran: beta = 0.95 ** 0.25 ≈ 0.9873
    beta: float = field(default_factory=lambda: 0.95 ** 0.25)
    
    # =====================
    # NON-CONVEX Adjustment Costs (CRITICAL!)
    # =====================
    # Capital adjustment costs
    capirrev: float = 0.339      # Capital irreversibility cost (partial)
    capfix: float = 0.0          # Fixed capital adjustment cost
    
    # Labor adjustment costs  
    hirelin: float = field(default_factory=lambda: 0.018 * 4.0)  # Hiring cost = 0.072
    firelin: float = field(default_factory=lambda: 0.018 * 4.0)  # Firing cost = 0.072
    labfix: float = field(default_factory=lambda: 0.024 * 4.0)   # Fixed labor adj cost = 0.096
    
    # =====================
    # Uncertainty Process
    # =====================
    ajump: float = 1.60569110682638       # Volatility jump multiplier (aggregate)
    zjump: float = 4.11699578856773       # Volatility jump multiplier (idiosyncratic)
    uncpers: float = 0.940556523390567    # Uncertainty shock persistence
    uncfreq: float = 0.0257892725462263   # Conditional probability of unc shock
    
    # =====================
    # Productivity Processes
    # =====================
    # Idiosyncratic productivity
    rhoz: float = 0.95
    sigmaz: float = 0.0507515557155377
    
    # Aggregate productivity  
    rhoa: float = 0.95
    sigmaa: float = 0.00668420914017636
    
    # =====================
    # Grid Sizes (simplified mode for faster computation)
    # =====================
    znum: int = 5          # Idiosyncratic productivity states (full: 9)
    anum: int = 7          # Aggregate productivity states (full: 21)
    snum: int = 2          # Volatility states (always 2)
    knum: int = 30         # Capital grid (full: 150)
    lnum: int = 15         # Labor grid (full: 75)
    kbarnum: int = 2       # Aggregate capital forecast grid
    
    # =====================
    # Grid Bounds
    # =====================
    kmin: float = 0.75
    lmin: float = 0.02
    amin: float = 0.8
    amax: float = 1.2
    kbarmin: float = 3.0
    kbarmax: float = 10.0
    
    # =====================
    # Simulation Controls
    # =====================
    Ncountries: int = 500      # Number of countries
    Tper: int = 100            # Time periods per country
    numdiscard: int = 500      # Burn-in periods
    nfirms: int = 800          # Total firms
    nfirmspub: int = 200       # Public firms (stock market)
    
    # =====================
    # VFI Controls
    # =====================
    vfmaxit: int = 50          # Max VFI iterations
    vferrortol: float = 1e-4   # VFI convergence tolerance
    accelmaxit: int = 200      # Howard acceleration iterations
    
    # =====================
    # Price/GE Controls
    # =====================
    maxpit: int = 50
    perrortol: float = 1e-3
    pval: float = 1.34         # Price level
    
    # =====================
    # IRF Controls
    # =====================
    numsimIRF: int = 2500      # Number of IRF simulations
    lengthIRF: int = 100       # Length of each IRF
    shockperIRF: int = 45      # Period to shock
    shocklengthIRF: int = 5    # Duration of low unc
    numdiscIRF: int = 45       # IRF burn-in
    
    # =====================
    # Disaster Parameters (for GMM)
    # =====================
    # Initial values from Fortran (best fit from uniform exploration)
    disaster_levels: Dict[str, float] = field(default_factory=lambda: {
        'nat_dis': -0.034,
        'pol_shock': 0.054,
        'revolution': -41.486,
        'terrorist': -3.950
    })
    
    disaster_unc_probs: Dict[str, float] = field(default_factory=lambda: {
        'nat_dis': 0.014,
        'pol_shock': 0.853,
        'revolution': 0.816,
        'terrorist': 0.110
    })
    
    # Parameter bounds for PSO
    param_bounds: Dict[str, Tuple[float, float]] = field(default_factory=lambda: {
        'nat_dis_levels': (-1.5, 0.1),
        'pol_shock_levels': (-0.1, 0.25),
        'revolution_levels': (-75.0, -35.5),
        'terrorist_levels': (-4.5, -2.5),
        'nat_dis_unc': (0.001, 0.35),
        'pol_shock_unc': (0.5, 0.999),
        'revolution_unc': (0.75, 0.999),
        'terrorist_unc': (0.001, 0.25)
    })
    
    # =====================
    # Data Moments (from paper)
    # =====================
    # First stage - macro sample
    first_stage_macro: np.ndarray = field(default_factory=lambda: np.array([
        [-0.071, -0.028],   # Nat Dis
        [1.657, 1.693],     # Pol Shock
        [-6.154, 7.841],    # Revolution
        [-0.047, -0.011]    # Terrorist
    ]))
    
    # Second stage - macro sample
    second_stage_macro: np.ndarray = field(default_factory=lambda: np.array([1.557, -3.859]))
    
    # First stage - micro sample
    first_stage_micro: np.ndarray = field(default_factory=lambda: np.array([
        [-0.147, 0.004],    # Nat Dis
        [1.852, 0.508],     # Pol Shock
        [-4.818, 3.201],    # Revolution
        [-0.117, 0.133]     # Terrorist
    ]))
    
    # Second stage - micro sample
    second_stage_micro: np.ndarray = field(default_factory=lambda: np.array([0.736, -9.735]))
    
    # Standard errors
    data_se: np.ndarray = field(default_factory=lambda: np.array([
        # First stage, levels LHS, macro
        0.1059, 0.0551, 1.0835, 0.0514,
        # First stage, vol LHS, macro
        0.0823, 0.1159, 2.2358, 0.0490,
        # Second stage, growth LHS, macro
        0.2906, 0.2844,
        # First stage, levels LHS, micro
        0.112, 0.085, 1.198, 0.044,
        # First stage, vol LHS, micro
        0.102, 0.130, 1.275, 0.083,
        # Second stage, growth LHS, micro
        0.558, 1.533
    ]))
    
    # Disaster probabilities (from data)
    disaster_probs: np.ndarray = field(default_factory=lambda: np.array([0.242, 0.03, 0.011, 0.008]))
    
    def __post_init__(self):
        """Compute derived quantities after initialization."""
        # Compute grid bounds based on parameters
        self.zmin = np.exp(-2.5 * np.sqrt(self.sigmaz**2 / (1 - self.rhoz**2)))
        self.zmax = np.exp(2.5 * np.sqrt(self.sigmaz**2 / (1 - self.rhoz**2)))
        
        # Compute high uncertainty ergodic probability
        self.highuncerg = self.uncfreq / (1.0 + self.uncfreq - self.uncpers)
    
    def get_data_moments(self) -> np.ndarray:
        """Return combined data moments vector (20 moments)."""
        moms = np.zeros(20)
        moms[0:4] = self.first_stage_macro[:, 0]    # First moment, macro
        moms[4:8] = self.first_stage_macro[:, 1]    # Second moment, macro
        moms[8:10] = self.second_stage_macro        # Second stage, macro
        moms[10:14] = self.first_stage_micro[:, 0]  # First moment, micro
        moms[14:18] = self.first_stage_micro[:, 1]  # Second moment, micro
        moms[18:20] = self.second_stage_micro       # Second stage, micro
        return moms
    
    def get_param_vector(self) -> np.ndarray:
        """Return disaster parameters as vector for optimization."""
        return np.array([
            self.disaster_levels['nat_dis'],
            self.disaster_levels['pol_shock'],
            self.disaster_levels['revolution'],
            self.disaster_levels['terrorist'],
            self.disaster_unc_probs['nat_dis'],
            self.disaster_unc_probs['pol_shock'],
            self.disaster_unc_probs['revolution'],
            self.disaster_unc_probs['terrorist']
        ])
    
    def set_param_vector(self, x: np.ndarray):
        """Set disaster parameters from optimization vector."""
        self.disaster_levels['nat_dis'] = x[0]
        self.disaster_levels['pol_shock'] = x[1]
        self.disaster_levels['revolution'] = x[2]
        self.disaster_levels['terrorist'] = x[3]
        self.disaster_unc_probs['nat_dis'] = x[4]
        self.disaster_unc_probs['pol_shock'] = x[5]
        self.disaster_unc_probs['revolution'] = x[6]
        self.disaster_unc_probs['terrorist'] = x[7]
    
    def get_param_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return lower and upper bounds for optimization."""
        lb = np.array([
            self.param_bounds['nat_dis_levels'][0],
            self.param_bounds['pol_shock_levels'][0],
            self.param_bounds['revolution_levels'][0],
            self.param_bounds['terrorist_levels'][0],
            self.param_bounds['nat_dis_unc'][0],
            self.param_bounds['pol_shock_unc'][0],
            self.param_bounds['revolution_unc'][0],
            self.param_bounds['terrorist_unc'][0]
        ])
        ub = np.array([
            self.param_bounds['nat_dis_levels'][1],
            self.param_bounds['pol_shock_levels'][1],
            self.param_bounds['revolution_levels'][1],
            self.param_bounds['terrorist_levels'][1],
            self.param_bounds['nat_dis_unc'][1],
            self.param_bounds['pol_shock_unc'][1],
            self.param_bounds['revolution_unc'][1],
            self.param_bounds['terrorist_unc'][1]
        ])
        return lb, ub


def create_params(simplified: bool = True) -> ModelParameters:
    """
    Create parameter set with appropriate grid sizes.
    
    Parameters
    ----------
    simplified : bool
        If True, use smaller grids for faster computation.
        If False, use full grid sizes matching Fortran.
    """
    if simplified:
        return ModelParameters()
    else:
        return ModelParameters(
            znum=9,
            anum=21,
            knum=150,
            lnum=75
        )
