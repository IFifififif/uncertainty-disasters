#!/usr/bin/env python3
"""
Lightweight cross-module consistency check.

This script runs a small, reproducible subset of each module and writes:
1) output/quick_check_metrics.json
2) output/quick_check_discrepancy.csv
3) output/quick_check_discrepancy.md
"""

from __future__ import annotations

import json
import time
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.iv.panel_iv import PanelIV
from src.iv_var.estimation import IVVAR
from src.lmn_var.estimation import LMNVAR
from src.model.solve import MicroMacroModel
from src.model.gmm import compute_simulated_moments

OUTPUT_DIR = PROJECT_ROOT / "output"
METRICS_JSON = OUTPUT_DIR / "quick_check_metrics.json"
DISCREPANCY_CSV = OUTPUT_DIR / "quick_check_discrepancy.csv"
DISCREPANCY_MD = OUTPUT_DIR / "quick_check_discrepancy.md"


def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return None


def run_iv() -> dict:
    t0 = time.time()
    iv = PanelIV()
    iv.load_data()
    t2 = iv.table2_baseline()
    col2 = t2["col2_iv"]
    return {
        "ok": True,
        "secs": round(time.time() - t0, 2),
        "coef_col2": [_safe_float(col2["coef_endog"][0]), _safe_float(col2["coef_endog"][1])],
        "se_col2": [_safe_float(col2["se_endog"][0]), _safe_float(col2["se_endog"][1])],
        "Jp_col2": _safe_float(col2.get("J_pval")),
        "F1_col2": _safe_float(col2["first_stage"][0]["F_stat"]),
    }


def run_iv_var() -> dict:
    t0 = time.time()
    ivv = IVVAR()
    ivv.load_data()
    baseline = ivv.estimate_baseline(seed=3991)
    # Small bootstrap for quick check speed; full replication still uses 150.
    se = ivv.bootstrap_se(baseline, n_boot=30, seed=3991, block_size=25)
    return {
        "ok": True,
        "secs": round(time.time() - t0, 2),
        "irf_first5": [float(x) for x in baseline["IRF_S_TO_Y"][:5]],
        "se_first5": [float(x) for x in se[:5]],
        "impact_t1": float(baseline["IRF_S_TO_Y"][0]),
    }


def run_lmn_var() -> dict:
    t0 = time.time()
    lmn = LMNVAR()
    lmn.load_data()
    lmn.step1_estimate_var_fe(include_country_fe=True, include_time_fe=True, nlags=12)
    out = lmn.step2_admissible_sets(
        n_draws=5000,
        seed=25041,
        restrictions={
            "rev3_min": 0.15,
            "rev2_max": -0.10,
            "coup3_min": 0.15,
            "coup2_min": 0.10,
            "nat2_max": 0.0,
            "ter2_max": 0.0,
        },
    )
    n_admit = int(out["n_admissible"])
    share = float(n_admit / 5000.0)
    return {
        "ok": True,
        "secs": round(time.time() - t0, 2),
        "n_admit": n_admit,
        "share": share,
        "impact_mean": float(out["IRF_med"][0, 2, 0]) if out.get("IRF_med") is not None else None,
        "impact_lb_t1": float(out["IRF_admit_lb"][0, 2, 0]) if out.get("IRF_admit_lb") is not None else None,
        "impact_ub_t1": float(out["IRF_admit_ub"][0, 2, 0]) if out.get("IRF_admit_ub") is not None else None,
        "impact_maxg_t1": float(out["IRF_maxg"][0, 2, 0]) if out.get("IRF_maxg") is not None else None,
    }


def run_model() -> dict:
    t0 = time.time()
    model = MicroMacroModel(simplified=True)
    model.build()
    model.solve(verbose=False)
    model.compute_irf(T=40, n_sims=100)

    p = model.params
    sim_mom = compute_simulated_moments(
        p,
        model.grids,
        model.vfi_solution,
        p.get_param_vector()[:4],
        p.get_param_vector()[4:],
    )
    data_mom = p.get_data_moments()
    rmse = float(np.sqrt(np.mean((sim_mom - data_mom) ** 2)))
    mae = float(np.mean(np.abs(sim_mom - data_mom)))

    return {
        "ok": True,
        "secs": round(time.time() - t0, 2),
        "irf_gdp_min": float(np.min(model.irf_results.irf_Y)),
        "irf_inv_min": float(np.min(model.irf_results.irf_I)),
        "mom_rmse": rmse,
        "mom_mae": mae,
    }


def build_discrepancy(metrics: dict) -> pd.DataFrame:
    rows = []
    iv = metrics.get("IV", {})
    ivv = metrics.get("IV_VAR", {})
    lmn = metrics.get("LMN_VAR", {})
    mdl = metrics.get("MODEL", {})

    if not iv.get("ok", False):
        rows.append({
            "module": "IV", "metric": "Execution status", "reference": "must run end-to-end",
            "python_value": 0.0, "delta": iv.get("error", "failed"), "status": "warning",
        })
        return pd.DataFrame(rows)

    rows.append({
        "module": "IV",
        "metric": "Table2 Col2 coef(cs_index_ret)",
        "reference": "Stata ivreg2 output",
        "python_value": iv.get("coef_col2", [None, None])[0],
        "delta": "N/A (no baseline file)",
        "status": "run-ok",
    })
    rows.append({
        "module": "IV",
        "metric": "Table2 Col2 coef(cs_index_vol)",
        "reference": "Stata ivreg2 output",
        "python_value": iv.get("coef_col2", [None, None])[1],
        "delta": "N/A (no baseline file)",
        "status": "run-ok",
    })
    rows.append({
        "module": "IV_VAR",
        "metric": "Baseline impact t=1",
        "reference": "MATLAB VAR.m (negative GDP response)",
        "python_value": ivv.get("impact_t1"),
        "delta": "sign-ok" if ivv.get("impact_t1") is not None and ivv["impact_t1"] < 0 else "sign-mismatch",
        "status": "run-ok" if ivv.get("impact_t1") is not None and ivv["impact_t1"] < 0 else "warning",
    })
    rows.append({
        "module": "IV_VAR",
        "metric": "Bootstrap SE finite check",
        "reference": "finite SE vector",
        "python_value": (ivv.get("se_first5") or [None, None, None, None])[3],
        "delta": "ok" if ivv.get("se_first5") and np.isfinite(ivv["se_first5"][3]) else "non-finite",
        "status": "run-ok" if ivv.get("se_first5") and np.isfinite(ivv["se_first5"][3]) else "warning",
    })
    rows.append({
        "module": "LMN_VAR",
        "metric": "Execution status",
        "reference": "must run end-to-end",
        "python_value": 1.0 if lmn.get("ok", False) else 0.0,
        "delta": "ok" if lmn.get("ok", False) else lmn.get("error", "failed"),
        "status": "run-ok" if lmn.get("ok", False) else "warning",
    })
    rows.append({
        "module": "LMN_VAR",
        "metric": "Admissible share (5k draws)",
        "reference": ">0 often requires large N draws",
        "python_value": lmn.get("share"),
        "delta": "ok" if lmn.get("share") is not None and lmn["share"] > 0 else "0 at low N draws",
        "status": "run-ok" if lmn.get("share") is not None and lmn["share"] > 0 else "warning",
    })
    rows.append({
        "module": "MODEL",
        "metric": "Moment RMSE (20 moments)",
        "reference": "0 (perfect match target)",
        "python_value": mdl.get("mom_rmse"),
        "delta": mdl.get("mom_rmse"),
        "status": "run-ok" if mdl.get("mom_rmse") is not None and mdl["mom_rmse"] < 10 else "warning",
    })
    rows.append({
        "module": "MODEL",
        "metric": "Moment MAE (20 moments)",
        "reference": "0 (perfect match target)",
        "python_value": mdl.get("mom_mae"),
        "delta": mdl.get("mom_mae"),
        "status": "run-ok" if mdl.get("mom_mae") is not None and mdl["mom_mae"] < 5 else "warning",
    })
    return pd.DataFrame(rows)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    metrics = {}
    for name, fn in [("IV", run_iv), ("IV_VAR", run_iv_var), ("LMN_VAR", run_lmn_var), ("MODEL", run_model)]:
        try:
            metrics[name] = fn()
        except Exception as exc:
            metrics[name] = {"ok": False, "error": str(exc)}

    with METRICS_JSON.open("w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2, ensure_ascii=False)

    disc = build_discrepancy(metrics)
    disc.to_csv(DISCREPANCY_CSV, index=False)

    with DISCREPANCY_MD.open("w", encoding="utf-8") as f:
        f.write("# Quick Check Discrepancy Table\n\n")
        f.write(disc.to_markdown(index=False))
        f.write("\n")

    print(f"Wrote {METRICS_JSON}")
    print(f"Wrote {DISCREPANCY_CSV}")
    print(f"Wrote {DISCREPANCY_MD}")


if __name__ == "__main__":
    main()
