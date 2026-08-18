#!/usr/bin/env python3
"""D1 minimum merge and first-stage diagnostics.

This script keeps the original BBT data untouched. It downloads public annual
buffer variables, maps BBT country codes to ISO3, merges lagged annual buffers
onto the BBT quarterly IV panel, and tests whether interacted disaster
instruments retain first-stage relevance.
"""

from __future__ import annotations

import io
import json
import math
import zipfile
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
import statsmodels.api as sm
from linearmodels.iv import IV2SLS
from scipy.linalg import qr


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "d1_validation"
RAW = OUT / "raw"

PANEL_PATH = ROOT / "data" / "IV" / "panel_iv_data.dta"
FISCAL_XLSX_URL = (
    "https://thedocs.worldbank.org/en/doc/"
    "31a9f7ffd7c72b9fbed93f3fb79af70a-0050012026/original/"
    "Fiscal-space-data.xlsx"
)

WB_INDICATORS = {
    "domestic_credit_private_gdp": "FS.AST.PRVT.GD.ZS",
    "gfdd_nonlife_insurance_prem_gdp": "GFDD.DI.10",
}

BBT_TO_ISO3 = {
    "AUS": "AUS",
    "AUT": "AUT",
    "BEL": "BEL",
    "BNG": "BGD",
    "BRA": "BRA",
    "CAN": "CAN",
    "CHI": "CHN",
    "CHL": "CHL",
    "COL": "COL",
    "CZH": "CZE",
    "DEN": "DNK",
    "ECU": "ECU",
    "EGY": "EGY",
    "FIN": "FIN",
    "FRA": "FRA",
    "GBR": "GBR",
    "GER": "DEU",
    "GRE": "GRC",
    "HUN": "HUN",
    "IND": "IND",
    "INO": "IDN",
    "IRE": "IRL",
    "IRN": "IRN",
    "ISR": "ISR",
    "ITA": "ITA",
    "JAP": "JPN",
    "KEN": "KEN",
    "KOR": "KOR",
    "KUW": "KWT",
    "LUX": "LUX",
    "MAL": "MYS",
    "MEX": "MEX",
    "MOR": "MAR",
    "NEH": "NLD",
    "NIG": "NGA",
    "NOR": "NOR",
    "NWZ": "NZL",
    "PAK": "PAK",
    "PER": "PER",
    "PHI": "PHL",
    "POL": "POL",
    "POR": "PRT",
    "ROM": "ROU",
    "RUS": "RUS",
    "SAF": "ZAF",
    "SAU": "SAU",
    "SER": "SRB",
    "SNG": "SGP",
    "SPA": "ESP",
    "SWE": "SWE",
    "SWI": "CHE",
    "TAI": "TWN",
    "THI": "THA",
    "TUN": "TUN",
    "TUR": "TUR",
    "UKR": "UKR",
    "USA": "USA",
    "VEN": "VEN",
    "VIE": "VNM",
}

INSTRUMENTS = [
    "l1savgnatshock",
    "l1savgpolshock",
    "l1savgrevshock",
    "l1savgtershock",
]
ENDOGENOUS = ["cs_index_ret", "cs_index_vol"]
BASE_VARS = ["ydgdp", "country", "year", "quarter", "yq_int"] + ENDOGENOUS + INSTRUMENTS


def fetch_url(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())
    return dest


def fetch_world_bank_indicator(indicator: str, name: str) -> pd.DataFrame:
    cache = RAW / f"wb_{indicator}.json"
    if not cache.exists():
        all_rows = []
        page = 1
        while True:
            url = (
                "https://api.worldbank.org/v2/country/all/indicator/"
                f"{quote(indicator)}?format=json&per_page=20000&page={page}"
            )
            req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urlopen(req, timeout=120) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            if not isinstance(data, list) or len(data) < 2:
                break
            meta, rows = data
            all_rows.extend(rows)
            if page >= int(meta.get("pages", page)):
                break
            page += 1
        cache.write_text(json.dumps(all_rows), encoding="utf-8")
    rows = json.loads(cache.read_text(encoding="utf-8"))
    out = []
    for r in rows:
        value = r.get("value")
        if value is None:
            continue
        country = r.get("countryiso3code")
        if not country:
            continue
        out.append(
            {
                "iso3": country,
                "year": int(r["date"]),
                name: float(value),
            }
        )
    return pd.DataFrame(out)


def read_fiscal_space() -> tuple[pd.DataFrame, dict]:
    xlsx = fetch_url(FISCAL_XLSX_URL, RAW / "Fiscal-space-data.xlsx")
    xl = pd.ExcelFile(xlsx)
    notes = {"workbook_sheets": xl.sheet_names, "selected_columns": {}}
    preferred_sheets = {
        "ggdy": "fiscal_gross_debt_gdp",
        "pby": "fiscal_primary_balance_gdp",
        "fby": "fiscal_balance_gdp",
    }
    frames = []
    for sheet, value_name in preferred_sheets.items():
        if sheet not in xl.sheet_names:
            continue
        df = pd.read_excel(xlsx, sheet_name=sheet)
        if "Country Code" not in df.columns:
            continue
        year_cols = []
        for c in df.columns:
            try:
                year = int(float(c))
            except (TypeError, ValueError):
                continue
            if 1900 <= year <= 2100:
                year_cols.append(c)
        if not year_cols:
            continue
        tmp = df[["Country Code"] + year_cols].copy()
        tmp = tmp.melt(id_vars=["Country Code"], var_name="year", value_name=value_name)
        tmp = tmp.rename(columns={"Country Code": "iso3"})
        tmp["iso3"] = tmp["iso3"].astype(str).str.upper().str.strip()
        tmp["year"] = tmp["year"].astype(float).astype(int)
        tmp[value_name] = pd.to_numeric(tmp[value_name], errors="coerce")
        tmp = tmp.dropna(subset=[value_name])
        notes["selected_columns"][sheet] = [value_name]
        frames.append(tmp)
    if frames:
        fiscal = frames[0]
        for f in frames[1:]:
            fiscal = fiscal.merge(f, on=["iso3", "year"], how="outer")
        return fiscal, notes

    # Fallback for future workbook layouts that may use long-form sheets.
    for sheet in xl.sheet_names:
        try:
            df = pd.read_excel(xlsx, sheet_name=sheet)
        except Exception:
            continue
        cols_lower = {str(c).strip().lower(): c for c in df.columns}
        country_col = None
        year_col = None
        for key, col in cols_lower.items():
            if key in {"code", "countrycode", "country code", "iso3", "iso"}:
                country_col = col
            if key == "year":
                year_col = col
        if country_col is None or year_col is None:
            continue
        candidate_cols = []
        for c in df.columns:
            cl = str(c).strip().lower()
            if any(token in cl for token in ["debt", "fiscal balance", "overall balance", "primary balance"]):
                candidate_cols.append(c)
        if not candidate_cols:
            continue
        keep = [country_col, year_col] + candidate_cols[:8]
        tmp = df[keep].copy()
        tmp = tmp.rename(columns={country_col: "iso3", year_col: "year"})
        tmp["iso3"] = tmp["iso3"].astype(str).str.upper().str.strip()
        tmp["year"] = pd.to_numeric(tmp["year"], errors="coerce")
        tmp = tmp.dropna(subset=["iso3", "year"])
        tmp["year"] = tmp["year"].astype(int)
        for c in candidate_cols[:8]:
            new = f"fiscal_{str(c).strip().lower().replace(' ', '_').replace('/', '_')}"
            tmp = tmp.rename(columns={c: new})
        notes["selected_columns"][sheet] = [c for c in tmp.columns if c not in {"iso3", "year"}]
        frames.append(tmp)
    if not frames:
        return pd.DataFrame(columns=["iso3", "year"]), notes
    fiscal = frames[0]
    for f in frames[1:]:
        fiscal = fiscal.merge(f, on=["iso3", "year"], how="outer")
    fiscal = fiscal.loc[:, ~fiscal.columns.duplicated()]
    fiscal = fiscal.replace(["..", "NA", "N/A", ""], np.nan)
    for c in fiscal.columns:
        if c not in {"iso3", "year"}:
            fiscal[c] = pd.to_numeric(fiscal[c], errors="coerce")
    return fiscal, notes


def standardize_within_sample(s: pd.Series) -> pd.Series:
    sd = s.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean()) / sd


def residualize(values: np.ndarray, design: pd.DataFrame) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    d = np.asarray(design, dtype=float)
    beta = np.linalg.lstsq(d, x, rcond=None)[0]
    return x - d @ beta


def cluster_f_test(y: pd.Series, x_cols: list[str], z_cols: list[str], data: pd.DataFrame) -> dict:
    d = data[[y.name] + x_cols + z_cols + ["country"]].dropna().copy()
    if d.empty:
        return {"nobs": 0, "nclusters": 0, "f_stat": math.nan, "p_value": math.nan}
    X = sm.add_constant(d[x_cols + z_cols], has_constant="add")
    model = sm.OLS(d[y.name], X).fit(cov_type="cluster", cov_kwds={"groups": d["country"]})
    terms = " = 0, ".join(z_cols) + " = 0"
    test = model.f_test(terms)
    return {
        "nobs": int(model.nobs),
        "nclusters": int(d["country"].nunique()),
        "f_stat": float(np.asarray(test.fvalue).ravel()[0]),
        "p_value": float(np.asarray(test.pvalue).ravel()[0]),
        "r2": float(model.rsquared),
    }


def first_stage_for_buffer(df: pd.DataFrame, buffer_col: str) -> list[dict]:
    needed = BASE_VARS + [buffer_col]
    d = df[needed].dropna().copy()
    if d.empty:
        return []
    d[f"{buffer_col}_z"] = standardize_within_sample(d[buffer_col])
    for z in INSTRUMENTS:
        d[f"{z}_x_{buffer_col}"] = d[z] * d[f"{buffer_col}_z"]
    for e in ENDOGENOUS:
        d[f"{e}_x_{buffer_col}"] = d[e] * d[f"{buffer_col}_z"]

    fe = pd.get_dummies(d["country"], prefix="c", drop_first=True, dtype=float)
    te = pd.get_dummies(d["yq_int"], prefix="t", drop_first=True, dtype=float)
    design = pd.concat([fe, te], axis=1)
    residualized = pd.DataFrame(index=d.index)
    regressors = ENDOGENOUS + [f"{e}_x_{buffer_col}" for e in ENDOGENOUS]
    instruments = INSTRUMENTS + [f"{z}_x_{buffer_col}" for z in INSTRUMENTS]
    for col in regressors + instruments:
        residualized[col] = residualize(d[col].to_numpy(), design)
    residualized["country"] = d["country"].values

    rows = []
    for target in regressors:
        base_test = cluster_f_test(
            residualized[target],
            x_cols=[],
            z_cols=INSTRUMENTS,
            data=residualized,
        )
        interaction_test = cluster_f_test(
            residualized[target],
            x_cols=INSTRUMENTS,
            z_cols=[f"{z}_x_{buffer_col}" for z in INSTRUMENTS],
            data=residualized,
        )
        all_test = cluster_f_test(
            residualized[target],
            x_cols=[],
            z_cols=instruments,
            data=residualized,
        )
        rows.append(
            {
                "buffer": buffer_col,
                "target_endogenous": target,
                "nobs": all_test["nobs"],
                "nclusters": all_test["nclusters"],
                "base_instr_F": base_test["f_stat"],
                "base_instr_p": base_test["p_value"],
                "interaction_instr_F": interaction_test["f_stat"],
                "interaction_instr_p": interaction_test["p_value"],
                "all_excluded_instr_F": all_test["f_stat"],
                "all_excluded_instr_p": all_test["p_value"],
                "all_excluded_instr_r2": all_test["r2"],
            }
        )
    return rows


def second_stage_for_buffer(df: pd.DataFrame, buffer_col: str) -> list[dict]:
    needed = BASE_VARS + [buffer_col]
    d = df[needed].dropna().copy()
    if d.empty:
        return []
    d[f"{buffer_col}_z"] = standardize_within_sample(d[buffer_col])
    for z in INSTRUMENTS:
        d[f"{z}_x_{buffer_col}"] = d[z] * d[f"{buffer_col}_z"]
    for e in ENDOGENOUS:
        d[f"{e}_x_{buffer_col}"] = d[e] * d[f"{buffer_col}_z"]

    fe = pd.get_dummies(d["country"], prefix="c", drop_first=True, dtype=float)
    te = pd.get_dummies(d["yq_int"], prefix="t", drop_first=True, dtype=float)
    design = pd.concat([fe, te], axis=1)
    endog_cols = ENDOGENOUS + [f"{e}_x_{buffer_col}" for e in ENDOGENOUS]
    instr_cols = INSTRUMENTS + [f"{z}_x_{buffer_col}" for z in INSTRUMENTS]
    resid = pd.DataFrame(index=d.index)
    resid["ydgdp"] = residualize(d["ydgdp"].to_numpy(), design)
    for col in endog_cols + instr_cols:
        resid[col] = residualize(d[col].to_numpy(), design)

    z = resid[instr_cols].copy()
    _, r, piv = qr(z.to_numpy(), mode="economic", pivoting=True)
    if r.size == 0:
        keep_instr = []
    else:
        tol = np.finfo(float).eps * max(z.shape) * abs(r[0, 0])
        rank = int((np.abs(np.diag(r)) > tol).sum())
        keep_instr = [z.columns[i] for i in sorted(piv[:rank])]
    if len(keep_instr) < len(endog_cols):
        raise ValueError(f"Only {len(keep_instr)} independent instruments for {len(endog_cols)} endogenous terms")

    model = IV2SLS(
        dependent=resid["ydgdp"],
        exog=None,
        endog=resid[endog_cols],
        instruments=resid[keep_instr],
    )
    res = model.fit(cov_type="clustered", clusters=d["country"])
    rows = []
    for term in endog_cols:
        rows.append(
            {
                "buffer": buffer_col,
                "term": term,
                "nobs": int(res.nobs),
                "nclusters": int(d["country"].nunique()),
                "coef": float(res.params[term]),
                "std_error": float(res.std_errors[term]),
                "p_value": float(res.pvalues[term]),
                "instruments_total": len(instr_cols),
                "instruments_kept": len(keep_instr),
            }
        )
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)

    panel = pd.read_stata(PANEL_PATH, convert_categoricals=False)
    panel["iso3"] = panel["country"].map(BBT_TO_ISO3)
    panel["year"] = panel["year"].astype(int)
    missing_iso = sorted(panel.loc[panel["iso3"].isna(), "country"].dropna().unique().tolist())

    buffers = []
    source_notes = {"world_bank_indicators": WB_INDICATORS, "fiscal_space_url": FISCAL_XLSX_URL}
    for name, indicator in WB_INDICATORS.items():
        buffers.append(fetch_world_bank_indicator(indicator, name))
    fiscal, fiscal_notes = read_fiscal_space()
    source_notes["fiscal_space_notes"] = fiscal_notes
    if not fiscal.empty:
        preferred = []
        for c in fiscal.columns:
            cl = c.lower()
            if c in {"iso3", "year"}:
                continue
            if "gross" in cl and "debt" in cl and len(preferred) < 1:
                preferred.append(c)
            elif "overall" in cl and "balance" in cl and len(preferred) < 2:
                preferred.append(c)
            elif "primary" in cl and "balance" in cl and len(preferred) < 3:
                preferred.append(c)
        if not preferred:
            nonkey = [c for c in fiscal.columns if c not in {"iso3", "year"}]
            preferred = nonkey[:3]
        fiscal = fiscal[["iso3", "year"] + preferred].copy()
        fiscal = fiscal.rename(columns={c: f"fiscal_space_{i+1}_{c}" for i, c in enumerate(preferred)})
        buffers.append(fiscal)

    annual = None
    for b in buffers:
        annual = b if annual is None else annual.merge(b, on=["iso3", "year"], how="outer")
    annual = annual.drop_duplicates(["iso3", "year"])
    annual.to_csv(OUT / "annual_buffers_raw.csv", index=False)

    annual_lag = annual.copy()
    annual_lag["year"] = annual_lag["year"] + 1
    rename = {c: f"l1_{c}" for c in annual_lag.columns if c not in {"iso3", "year"}}
    annual_lag = annual_lag.rename(columns=rename)
    merged = panel.merge(annual_lag, on=["iso3", "year"], how="left")
    merged.to_csv(OUT / "bbt_d1_merged_panel_preview.csv", index=False)

    buffer_cols = list(rename.values())
    base_sample = panel[BASE_VARS].dropna()
    coverage_rows = []
    for col in buffer_cols:
        d = merged[BASE_VARS + [col]].dropna()
        coverage_rows.append(
            {
                "buffer": col,
                "nonmissing_panel_obs": int(merged[col].notna().sum()),
                "analysis_obs": int(len(d)),
                "analysis_countries": int(d["country"].nunique()) if len(d) else 0,
                "analysis_year_min": int(d["year"].min()) if len(d) else None,
                "analysis_year_max": int(d["year"].max()) if len(d) else None,
                "share_of_baseline_iv_sample": len(d) / len(base_sample) if len(base_sample) else math.nan,
            }
        )
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(OUT / "coverage_by_buffer.csv", index=False)

    fs_rows = []
    for col in buffer_cols:
        fs_rows.extend(first_stage_for_buffer(merged, col))
    first_stage = pd.DataFrame(fs_rows)
    first_stage.to_csv(OUT / "first_stage_diagnostics.csv", index=False)

    ss_rows = []
    for col in buffer_cols:
        try:
            ss_rows.extend(second_stage_for_buffer(merged, col))
        except Exception as exc:
            ss_rows.append(
                {
                    "buffer": col,
                    "term": "__ERROR__",
                    "nobs": 0,
                    "nclusters": 0,
                    "coef": math.nan,
                    "std_error": math.nan,
                    "p_value": math.nan,
                    "error": repr(exc),
                }
            )
    second_stage = pd.DataFrame(ss_rows)
    second_stage.to_csv(OUT / "second_stage_smoke_iv.csv", index=False)

    source_notes["missing_iso3_mappings"] = missing_iso
    source_notes["panel_rows"] = int(len(panel))
    source_notes["baseline_iv_complete_rows"] = int(len(base_sample))
    source_notes["output_files"] = [
        "annual_buffers_raw.csv",
        "bbt_d1_merged_panel_preview.csv",
        "coverage_by_buffer.csv",
        "first_stage_diagnostics.csv",
        "second_stage_smoke_iv.csv",
        "summary.json",
    ]
    if not first_stage.empty:
        source_notes["first_stage_summary"] = {
            "buffers_tested": sorted(first_stage["buffer"].unique().tolist()),
            "weak_interaction_tests_F_lt_10": int((first_stage["interaction_instr_F"] < 10).sum()),
            "strongest_interaction_F": float(first_stage["interaction_instr_F"].max()),
            "weakest_interaction_F": float(first_stage["interaction_instr_F"].min()),
        }
    (OUT / "summary.json").write_text(json.dumps(source_notes, indent=2, ensure_ascii=False), encoding="utf-8")

    print("Wrote D1 validation outputs to", OUT)
    print("\nCoverage:")
    print(coverage.to_string(index=False))
    print("\nFirst-stage diagnostics:")
    if first_stage.empty:
        print("No first-stage diagnostics produced.")
    else:
        print(first_stage.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\nSecond-stage IV smoke test:")
    if second_stage.empty:
        print("No second-stage diagnostics produced.")
    else:
        print(second_stage.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
