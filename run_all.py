"""
Main runner script for BBT (2024) replication.

Run all modules or individual modules:
    python run_all.py              # Run everything
    python run_all.py iv           # Tables 1-6
    python run_all.py iv_var       # Figures 6-7
    python run_all.py lmn_var      # Figures 3-5
    python run_all.py model        # Figure 8
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


def run_iv():
    """Run Panel IV regressions (Tables 1-6)."""
    from src.iv.panel_iv import PanelIV
    panel_iv = PanelIV()
    panel_iv.run_all()


def run_iv_var():
    """Run IV-VAR estimation (Figures 6-7)."""
    from src.iv_var.estimation import IVVAR
    ivvar = IVVAR()
    ivvar.run_all()


def run_lmn_var():
    """Run LMN VAR estimation (Figures 3-5)."""
    from src.lmn_var.estimation import LMNVAR
    lmn = LMNVAR()
    lmn.run_all()


def run_model():
    """Run model simulation (Figure 8)."""
    from src.model.solve import MicroMacroModel
    model = MicroMacroModel()
    model.run_all()


def main():
    modules = sys.argv[1:] if len(sys.argv) > 1 else ['iv', 'iv_var', 'lmn_var', 'model']

    runners = {
        'iv': run_iv,
        'iv_var': run_iv_var,
        'lmn_var': run_lmn_var,
        'model': run_model,
    }

    for mod in modules:
        if mod in runners:
            print(f"\n{'=' * 70}")
            print(f"Running module: {mod}")
            print(f"{'=' * 70}")
            runners[mod]()
        else:
            print(f"Unknown module: {mod}")
            print(f"Available: {list(runners.keys())}")
            sys.exit(1)

    print("\n" + "=" * 70)
    print("ALL MODULES COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
