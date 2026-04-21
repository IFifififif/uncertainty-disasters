#!/usr/bin/env python3
"""Run MODEL-vs-paper comparison with checkpoint-friendly options.

Examples
--------
python tools/model_paper_compare.py --profile medium --n-sims 20 --save-vfi output/model_cache/medium_vfi.npz
python tools/model_paper_compare.py --profile medium --n-sims 20 --load-vfi output/model_cache/medium_vfi.npz
python tools/model_paper_compare.py --profile full --n-sims 5 --irf-t 40
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import sys
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.model.solve import MicroMacroModel
from src.model.vfi import VFISolution


def apply_profile(model: MicroMacroModel, profile: str) -> None:
    """Apply reproducible grid profiles."""
    if profile == "full":
        # Fortran baseline grid (already default when simplified=False).
        return
    if profile == "medium":
        # Fidelity-oriented but runnable on laptop/CI.
        p = model.params
        p.znum = 7
        p.anum = 11
        p.knum = 40
        p.lnum = 20
        return
    if profile == "smoke":
        # Fast sanity profile, aligned with quick-check scale.
        p = model.params
        p.znum = 5
        p.anum = 7
        p.knum = 30
        p.lnum = 15
        return
    raise ValueError(f"Unknown profile: {profile}")


def apply_overrides(model: MicroMacroModel, znum: int | None, anum: int | None, knum: int | None, lnum: int | None) -> None:
    p = model.params
    if znum is not None:
        p.znum = int(znum)
    if anum is not None:
        p.anum = int(anum)
    if knum is not None:
        p.knum = int(knum)
    if lnum is not None:
        p.lnum = int(lnum)


def save_vfi(path: Path, sol: VFISolution, meta: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        V=sol.V,
        polmat=sol.polmat,
        kprime_pos=sol.kprime_pos,
        lpol_pos=sol.lpol_pos,
        EVmat=sol.EVmat,
        converged=np.array([int(sol.converged)], dtype=np.int32),
        iterations=np.array([int(sol.iterations)], dtype=np.int32),
        vf_error=np.array([float(sol.vf_error)], dtype=np.float64),
        pol_error=np.array([float(sol.pol_error)], dtype=np.float64),
        meta_json=np.array([json.dumps(meta, ensure_ascii=False)]),
    )


def load_vfi(path: Path) -> tuple[VFISolution, dict]:
    data = np.load(path, allow_pickle=True)
    meta = json.loads(str(data["meta_json"][0])) if "meta_json" in data else {}
    sol = VFISolution(
        V=data["V"],
        polmat=data["polmat"],
        kprime_pos=data["kprime_pos"],
        lpol_pos=data["lpol_pos"],
        EVmat=data["EVmat"],
        converged=bool(int(data["converged"][0])),
        iterations=int(data["iterations"][0]),
        vf_error=float(data["vf_error"][0]),
        pol_error=float(data["pol_error"][0]),
    )
    return sol, meta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=["full", "medium", "smoke"], default="medium")
    parser.add_argument("--irf-t", type=int, default=40)
    parser.add_argument("--n-sims", type=int, default=20)
    parser.add_argument("--seed", type=int, default=2501)
    parser.add_argument("--save-vfi", type=Path)
    parser.add_argument("--load-vfi", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--znum", type=int)
    parser.add_argument("--anum", type=int)
    parser.add_argument("--knum", type=int)
    parser.add_argument("--lnum", type=int)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    out = args.output or (PROJECT_ROOT / "output" / f"model_paper_compare_{args.profile}.json")
    out.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()

    # Keep solver path consistent across profiles: full VFI routine, profile only
    # changes grid sizes.
    model = MicroMacroModel(simplified=False)
    apply_profile(model, args.profile)
    apply_overrides(model, args.znum, args.anum, args.knum, args.lnum)

    t_build0 = time.time()
    model.build()
    t_build = time.time() - t_build0

    loaded_vfi = False
    load_meta = {}
    if args.load_vfi:
        model.vfi_solution, load_meta = load_vfi(args.load_vfi)
        loaded_vfi = True
    else:
        t_solve0 = time.time()
        model.solve(verbose=args.verbose)
        t_solve = time.time() - t_solve0

        if args.save_vfi:
            meta = {
                "profile": args.profile,
                "grid": {
                    "znum": model.params.znum,
                    "anum": model.params.anum,
                    "snum": model.params.snum,
                    "knum": model.params.knum,
                    "lnum": model.params.lnum,
                },
                "vf": {
                    "iterations": int(model.vfi_solution.iterations),
                    "converged": bool(model.vfi_solution.converged),
                },
                "saved_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            }
            save_vfi(args.save_vfi, model.vfi_solution, meta)
    if loaded_vfi:
        t_solve = None

    t_irf0 = time.time()
    model.compute_irf(T=args.irf_t, n_sims=args.n_sims)
    t_irf = time.time() - t_irf0

    irf_y = model.irf_results.irf_Y
    irf_i = model.irf_results.irf_I

    result = {
        "profile": args.profile,
        "seed": args.seed,
        "irf_t": args.irf_t,
        "n_sims": args.n_sims,
        "timing_sec": {
            "build": round(t_build, 3),
            "solve": None if t_solve is None else round(t_solve, 3),
            "irf": round(t_irf, 3),
            "total": round(time.time() - t0, 3),
        },
        "grid": {
            "znum": int(model.params.znum),
            "anum": int(model.params.anum),
            "snum": int(model.params.snum),
            "knum": int(model.params.knum),
            "lnum": int(model.params.lnum),
            "states": int(model.grids.numstates),
        },
        "vfi": {
            "loaded": loaded_vfi,
            "load_meta": load_meta,
            "converged": bool(model.vfi_solution.converged),
            "iterations": int(model.vfi_solution.iterations),
            "vf_error": float(model.vfi_solution.vf_error),
            "pol_error": float(model.vfi_solution.pol_error),
        },
        "irf": {
            "gdp_t1": float(irf_y[0]),
            "gdp_min": float(np.min(irf_y)),
            "inv_t1": float(irf_i[0]),
            "inv_min": float(np.min(irf_i)),
            "gdp_series": [float(x) for x in irf_y[: min(40, len(irf_y))]],
            "inv_series": [float(x) for x in irf_i[: min(40, len(irf_i))]],
        },
        # MODEL block is supplemental in the replication packet README; no exact
        # published paper target moments are provided for Figure 8 amplitudes.
        "paper_targets": {
            "gdp_drop_pp": None,
            "inv_drop_pp_approx": None,
        },
        "gap_vs_paper": {
            "gdp_min_gap": None,
            "inv_min_gap": None,
        },
    }

    with out.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    print(f"Wrote {out}")
    print(json.dumps({
        "profile": result["profile"],
        "states": result["grid"]["states"],
        "timing_sec": result["timing_sec"],
        "gdp_min": result["irf"]["gdp_min"],
        "inv_min": result["irf"]["inv_min"],
        "gap_vs_paper": result["gap_vs_paper"],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
