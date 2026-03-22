"""
Basic smoke tests for BBT replication project.

These tests verify that:
1. Data loads correctly
2. Modules can be instantiated
3. Basic estimation runs without errors
"""

import sys
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def test_data_loading():
    """Test that all data files load correctly."""
    import pandas as pd

    # IV data
    df_iv = pd.read_stata(PROJECT_ROOT / "data" / "IV" / "panel_iv_data.dta")
    assert len(df_iv) > 0, "IV data is empty"
    assert 'ydgdp' in df_iv.columns, "Missing ydgdp in IV data"
    print(f"  [PASS] IV data: {len(df_iv)} obs, {len(df_iv.columns)} cols")

    # IV_VAR data
    df_ivvar = pd.read_csv(PROJECT_ROOT / "data" / "IV_VAR" / "VARdata.csv")
    assert len(df_ivvar) > 0, "IV_VAR data is empty"
    assert df_ivvar.shape[1] >= 7, "IV_VAR data should have >= 7 columns"
    print(f"  [PASS] IV_VAR data: {df_ivvar.shape}")

    # LMN_VAR data
    df_lmn = pd.read_stata(PROJECT_ROOT / "data" / "LMN_VAR" / "Dates_and_Data.dta")
    assert len(df_lmn) > 0, "LMN_VAR data is empty"
    assert 'ydgdp' in df_lmn.columns, "Missing ydgdp in LMN_VAR data"
    print(f"  [PASS] LMN_VAR data: {len(df_lmn)} obs, {len(df_lmn.columns)} cols")


def test_iv_module():
    """Test Panel IV module instantiation and basic run."""
    from src.iv.panel_iv import PanelIV
    panel_iv = PanelIV()
    panel_iv.load_data()
    assert panel_iv.df is not None
    assert len(panel_iv.df) > 0
    print("  [PASS] IV module loads and instantiates")


def test_iv_var_module():
    """Test IV-VAR module instantiation."""
    from src.iv_var.estimation import IVVAR
    ivvar = IVVAR()
    ivvar.load_data()
    assert ivvar.data is not None
    print("  [PASS] IV-VAR module loads and instantiates")


def test_lmn_var_module():
    """Test LMN VAR module instantiation."""
    from src.lmn_var.estimation import LMNVAR
    lmn = LMNVAR()
    lmn.load_data()
    assert lmn.df is not None
    print("  [PASS] LMN VAR module loads and instantiates")


def test_model_module():
    """Test Model module instantiation."""
    from src.model.solve import MicroMacroModel
    model = MicroMacroModel()
    assert model.params is not None
    print("  [PASS] Model module instantiates")


def test_utils():
    """Test utility functions."""
    from src.utils.regression import ols_with_cluster_se, iv2sls_with_cluster_se

    # Simple OLS test
    np.random.seed(42)
    N = 100
    X = np.random.randn(N, 2)
    y = 1.5 * X[:, 0] + 0.5 * X[:, 1] + 0.1 * np.random.randn(N)
    clusters = np.repeat(np.arange(10), 10)

    result = ols_with_cluster_se(y, X, clusters)
    assert abs(result['coef'][0] - 1.5) < 0.3, f"OLS coef should be ~1.5, got {result['coef'][0]}"
    print(f"  [PASS] OLS: coef={result['coef']}, se={result['se']}")

    # Simple IV test
    Z = np.random.randn(N, 2)
    X_endog = X[:, :1]
    result_iv = iv2sls_with_cluster_se(y, X_endog, X[:, 1:], Z, clusters)
    assert result_iv['nobs'] == N
    print(f"  [PASS] IV2SLS: nobs={result_iv['nobs']}")


if __name__ == '__main__':
    print("=" * 50)
    print("Running smoke tests")
    print("=" * 50)

    test_data_loading()
    test_utils()
    test_iv_module()
    test_iv_var_module()
    test_lmn_var_module()
    test_model_module()

    print("\n" + "=" * 50)
    print("ALL TESTS PASSED")
    print("=" * 50)
