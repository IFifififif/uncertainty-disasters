"""
Unified runner for BBT replication with centralized JSON config.

Usage examples:
    python run_all.py
    python run_all.py --config config/experiment_config.json
    python run_all.py --jobs baseline alt_data
    python run_all.py iv iv_var --sequential
"""

import argparse
import concurrent.futures as cf
import json
import traceback
from pathlib import Path
from typing import Dict, List

PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "experiment_config.json"
MODULES = ["iv", "iv_var", "lmn_var", "model"]


def _resolve_path(path_str: str) -> str:
    p = Path(path_str)
    if p.is_absolute():
        return str(p)
    return str((PROJECT_ROOT / p).resolve())


def _load_config(config_path: Path) -> dict:
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _as_bool(value, key: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes", "y", "on"}:
            return True
        if lowered in {"false", "0", "no", "n", "off"}:
            return False
    raise ValueError(f"Config value `{key}` must be a boolean, got {value!r}")


def _merge_module_cfg(defaults: dict, job: dict, module: str) -> dict:
    # Keep documentation keys (for human readability in config files) harmless:
    # execution code only reads known runtime keys via .get(...).
    cfg = dict(defaults.get(module, {}))
    if isinstance(job.get(module), dict):
        cfg.update(job[module])
    return cfg


def _normalize_jobs(cfg: dict, selected_jobs: List[str]) -> List[dict]:
    jobs = cfg.get("jobs", [])
    if not jobs:
        raise ValueError("Config must include at least one job in `jobs`.")
    out = []
    selected = set(selected_jobs or [])
    for job in jobs:
        name = job.get("name")
        if not name:
            continue
        if not job.get("enabled", True):
            continue
        if selected and name not in selected:
            continue
        out.append(job)
    if not out:
        raise ValueError("No enabled jobs matched selection.")
    return out


def _run_one_job(payload: dict) -> dict:
    """
    Run one job (optionally in subprocess).
    Returns a compact status dict.
    """
    job_name = payload["job_name"]
    modules = payload["modules"]
    module_cfgs = payload["module_cfgs"]
    run_root = Path(payload["run_root"])
    run_root.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 78)
    print(f"Running job: {job_name}")
    print(f"Run directory: {run_root}")
    print(f"Modules: {modules}")
    print("=" * 78)

    for mod in modules:
        print(f"\n--- [{job_name}] module={mod} ---")
        if mod == "iv":
            from src.iv.panel_iv import PanelIV
            c = module_cfgs["iv"]
            iv = PanelIV(
                data_path=_resolve_path(c["data_path"]),
                standardize_residualized=c.get("standardize_residualized", "none"),
                output_dir=str(run_root / "tables"),
            )
            iv.run_all()

        elif mod == "iv_var":
            from src.iv_var.estimation import IVVAR
            c = module_cfgs["iv_var"]
            ivv = IVVAR(
                data_path=_resolve_path(c["data_path"]),
                output_dir=str(run_root / "figures"),
            )
            ivv.run_all(
                n_starts=int(c.get("n_starts", 1)),
                start_jitter=float(c.get("start_jitter", 0.0)),
                selection_mode=str(c.get("selection_mode", "objective")),
                objective_tie_tol=float(c.get("objective_tie_tol", 2e-12)),
                target_impact_t1=float(c.get("target_impact_t1", -3.5)),
                diag_floor=float(c.get("diag_floor", 0.16)),
                bootstrap_n=int(c.get("bootstrap_n", 150)),
                bootstrap_block_size=int(c.get("bootstrap_block_size", 25)),
            )

        elif mod == "lmn_var":
            from src.lmn_var.estimation import LMNVAR
            c = module_cfgs["lmn_var"]
            lmn = LMNVAR(
                data_path=_resolve_path(c["data_path"]),
                output_dir=str(run_root / "figures"),
            )
            lmn.run_all(n_draws=int(c.get("n_draws", 1500000)))

        elif mod == "model":
            from src.model.solve import MicroMacroModel
            c = module_cfgs["model"]
            model = MicroMacroModel(
                simplified=_as_bool(c.get("simplified", False), "model.simplified"),
                output_dir=str(run_root / "figures"),
            )
            model.run_all(
                do_estimation=_as_bool(c.get("do_estimation", False), "model.do_estimation"),
                irf_T=int(c.get("irf_t", 40)),
                irf_n_sims=int(c.get("irf_n_sims", 100)),
            )
        else:
            raise ValueError(f"Unknown module: {mod}")

    return {"job": job_name, "ok": True, "run_root": str(run_root)}


def _build_payloads(cfg: dict, jobs: List[dict], cli_modules: List[str]) -> List[dict]:
    defaults = cfg.get("defaults", {})
    cfg_default_modules = cfg.get("default_modules", MODULES)
    payloads: List[dict] = []
    for job in jobs:
        job_name = job["name"]
        modules = cli_modules if cli_modules else job.get("modules", cfg_default_modules)
        for m in modules:
            if m not in MODULES:
                raise ValueError(f"Unsupported module: {m}")

        module_cfgs = {m: _merge_module_cfg(defaults, job, m) for m in MODULES}
        run_root = PROJECT_ROOT / "output" / "runs" / job_name

        payloads.append(
            {
                "job_name": job_name,
                "modules": modules,
                "module_cfgs": module_cfgs,
                "run_root": str(run_root),
            }
        )
    return payloads


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("modules", nargs="*", help="Optional modules to run: iv iv_var lmn_var model")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH), help="Path to unified JSON config")
    parser.add_argument("--jobs", nargs="*", default=[], help="Optional job names to run")
    parser.add_argument("--sequential", action="store_true", help="Force sequential execution")
    args = parser.parse_args()

    cfg = _load_config(Path(args.config))
    jobs = _normalize_jobs(cfg, args.jobs)
    payloads = _build_payloads(cfg, jobs, args.modules)

    par_cfg = cfg.get("parallel", {})
    parallel_enabled = bool(par_cfg.get("enabled", False)) and not args.sequential and len(payloads) > 1
    max_workers = int(par_cfg.get("max_workers", 2))

    results = []
    if parallel_enabled:
        print(f"Running jobs in parallel with max_workers={max_workers}")
        with cf.ProcessPoolExecutor(max_workers=max_workers) as pool:
            fut_to_job = {pool.submit(_run_one_job, p): p["job_name"] for p in payloads}
            for fut in cf.as_completed(fut_to_job):
                job_name = fut_to_job[fut]
                try:
                    res = fut.result()
                    results.append(res)
                    print(f"[DONE] {job_name}")
                except Exception as e:
                    tb = traceback.format_exc()
                    results.append({"job": job_name, "ok": False, "error": str(e), "traceback": tb})
                    print(f"[FAILED] {job_name}: {e}")
    else:
        for p in payloads:
            job_name = p["job_name"]
            try:
                res = _run_one_job(p)
                results.append(res)
                print(f"[DONE] {job_name}")
            except Exception as e:
                tb = traceback.format_exc()
                results.append({"job": job_name, "ok": False, "error": str(e), "traceback": tb})
                print(f"[FAILED] {job_name}: {e}")

    summary_path = PROJECT_ROOT / "output" / "runs" / "run_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    n_fail = sum(0 if r.get("ok") else 1 for r in results)
    print("\n" + "=" * 78)
    print(f"Run complete. jobs={len(results)}, failed={n_fail}")
    print(f"Summary: {summary_path}")
    print("=" * 78)
    if n_fail > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
