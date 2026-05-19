#!/usr/bin/env python3
"""
Build a configurable LaTeX comparison report from run outputs.

The report content is intentionally driven by the JSON config instead of hard
coded paper-specific labels. Users can rename variables, choose which generated
CSV tables and PDF figures to include, and define numeric comparison checks.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def latex_escape(value: Any) -> str:
    text = "" if value is None else str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def resolve_path(path_value: str | None, base: Path = PROJECT_ROOT) -> Path | None:
    if not path_value:
        return None
    path = Path(path_value)
    return path if path.is_absolute() else (base / path).resolve()


def fmt_number(value: Any, digits: int = 6) -> str:
    if value is None or value == "":
        return ""
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def get_nested(data: Any, dotted_path: str) -> Any:
    cur = data
    for part in dotted_path.split("."):
        if isinstance(cur, dict):
            cur = cur[part]
        elif isinstance(cur, list):
            cur = cur[int(part)]
        else:
            raise KeyError(dotted_path)
    return cur


def row_matches(row: dict[str, str], filters: dict[str, Any]) -> bool:
    return all(str(row.get(k, "")) == str(v) for k, v in filters.items())


def load_metric_value(spec: dict[str, Any], output_root: Path) -> Any:
    source_type = spec.get("source_type", "csv")
    path = resolve_path(spec.get("path"), PROJECT_ROOT)
    if path is None:
        run_name = spec.get("run")
        rel_file = spec.get("file")
        if not run_name or not rel_file:
            raise ValueError(f"Comparison `{spec.get('name', '')}` needs path or run+file.")
        path = output_root / "runs" / run_name / rel_file

    if source_type == "json":
        data = load_json(path)
        return get_nested(data, spec["json_path"])

    if source_type == "csv":
        rows = read_csv_rows(path)
        filters = spec.get("row_filter", {})
        matches = [r for r in rows if row_matches(r, filters)]
        if not matches:
            raise ValueError(f"No CSV row matched {filters} in {path}")
        return matches[0][spec["value_column"]]

    raise ValueError(f"Unsupported source_type: {source_type}")


def assess(value: Any, expected: Any, tolerance: float | None, comparator: str) -> tuple[str, str]:
    if expected is None or expected == "":
        return "reported", ""
    try:
        actual = float(value)
        target = float(expected)
    except (TypeError, ValueError):
        return ("pass" if str(value) == str(expected) else "check"), ""

    delta = actual - target
    if comparator == "abs_tol":
        ok = tolerance is not None and abs(delta) <= tolerance
    elif comparator == "lt":
        ok = actual < target
    elif comparator == "gt":
        ok = actual > target
    elif comparator == "lte":
        ok = actual <= target
    elif comparator == "gte":
        ok = actual >= target
    else:
        raise ValueError(f"Unsupported comparator: {comparator}")
    return ("pass" if ok else "check"), fmt_number(delta)


def table_from_csv(spec: dict[str, Any], output_root: Path, labels: dict[str, str]) -> list[str]:
    run_name = spec.get("run")
    path = resolve_path(spec.get("path"), PROJECT_ROOT)
    if path is None:
        path = output_root / "runs" / run_name / spec["file"]
    title = spec.get("title", path.name)
    max_rows = int(spec.get("max_rows", 40))

    out = [r"\subsection*{" + latex_escape(title) + "}"]
    if not path.exists():
        out.append(latex_escape(f"Missing table output: {path}"))
        return out

    rows = read_csv_rows(path)
    columns = spec.get("columns") or (list(rows[0].keys()) if rows else [])
    if not rows:
        out.append(latex_escape(f"No rows found at {path}"))
        return out

    col_def = "l" * len(columns)
    out.append(r"\begin{longtable}{" + col_def + "}")
    out.append(r"\toprule")
    out.append(" & ".join(latex_escape(labels.get(c, c)) for c in columns) + r" \\")
    out.append(r"\midrule")
    out.append(r"\endfirsthead")
    out.append(r"\toprule")
    out.append(" & ".join(latex_escape(labels.get(c, c)) for c in columns) + r" \\")
    out.append(r"\midrule")
    out.append(r"\endhead")
    for row in rows[:max_rows]:
        out.append(" & ".join(latex_escape(row.get(c, "")) for c in columns) + r" \\")
    out.append(r"\bottomrule")
    out.append(r"\end{longtable}")
    if len(rows) > max_rows:
        out.append(latex_escape(f"Showing first {max_rows} of {len(rows)} rows."))
    return out


def build_tex(config: dict[str, Any], output_path: Path) -> None:
    report = config.get("report", {})
    output_root = resolve_path(report.get("output_root", "output"), PROJECT_ROOT) or (PROJECT_ROOT / "output")
    labels = report.get("labels", {})
    variable_labels = labels.get("variables", {})
    module_labels = labels.get("modules", {})
    date_value = report.get("date", r"\today")
    date_latex = date_value if date_value == r"\today" else latex_escape(date_value)

    lines: list[str] = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[a4paper,margin=1in]{geometry}",
        r"\usepackage{booktabs}",
        r"\usepackage{longtable}",
        r"\usepackage{array}",
        r"\usepackage{graphicx}",
        r"\usepackage{hyperref}",
        r"\usepackage{float}",
        r"\title{" + latex_escape(report.get("title", "Replication Comparison Report")) + "}",
        r"\author{" + latex_escape(report.get("author", "Local Python Pipeline")) + "}",
        r"\date{" + date_latex + "}",
        r"\begin{document}",
        r"\maketitle",
    ]

    if report.get("summary"):
        lines.extend([r"\section*{" + latex_escape(report.get("summary_title", "Summary")) + "}", latex_escape(report["summary"])])

    run_summary = resolve_path(report.get("run_summary", "output/runs/run_summary.json"), PROJECT_ROOT)
    if run_summary and run_summary.exists():
        lines.append(r"\section*{Run Status}")
        rows = load_json(run_summary)
        lines.extend([
            r"\begin{longtable}{llll}",
            r"\toprule",
            r"Job & Status & Run directory & Error \\",
            r"\midrule",
        ])
        for row in rows:
            status = "ok" if row.get("ok") else "failed"
            lines.append(
                " & ".join(
                    latex_escape(x)
                    for x in [
                        row.get("job", ""),
                        status,
                        row.get("run_root", ""),
                        row.get("error", ""),
                    ]
                )
                + r" \\"
            )
        lines.extend([r"\bottomrule", r"\end{longtable}"])

    comparisons = report.get("comparisons", [])
    if comparisons:
        lines.append(r"\section*{" + latex_escape(report.get("comparison_title", "Numeric Comparison")) + "}")
        lines.extend([
            r"\begin{longtable}{lllllll}",
            r"\toprule",
            r"Module & Metric & Expected & Actual & Delta & Tolerance & Status \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"Module & Metric & Expected & Actual & Delta & Tolerance & Status \\",
            r"\midrule",
            r"\endhead",
        ])
        for spec in comparisons:
            try:
                value = load_metric_value(spec, output_root)
                status, delta = assess(
                    value,
                    spec.get("expected"),
                    spec.get("tolerance"),
                    spec.get("comparator", "abs_tol"),
                )
                error = ""
            except Exception as exc:
                value = ""
                status = "error"
                delta = ""
                error = str(exc)
            module = module_labels.get(spec.get("module", ""), spec.get("module", ""))
            metric = spec.get("metric", spec.get("name", ""))
            if error:
                metric = f"{metric}: {error}"
            lines.append(
                " & ".join(
                    latex_escape(x)
                    for x in [
                        module,
                        metric,
                        fmt_number(spec.get("expected")),
                        fmt_number(value),
                        delta,
                        "" if spec.get("tolerance") is None else fmt_number(spec.get("tolerance")),
                        status,
                    ]
                )
                + r" \\"
            )
        lines.extend([r"\bottomrule", r"\end{longtable}"])

    tables = report.get("tables", [])
    if tables:
        lines.append(r"\section*{" + latex_escape(report.get("tables_title", "Generated Tables")) + "}")
        for spec in tables:
            lines.extend(table_from_csv(spec, output_root, variable_labels))

    figures = report.get("figures", [])
    if figures:
        lines.append(r"\section*{" + latex_escape(report.get("figures_title", "Generated Figures")) + "}")
        for spec in figures:
            path = resolve_path(spec.get("path"), PROJECT_ROOT)
            if path is None:
                path = output_root / "runs" / spec["run"] / spec["file"]
            title = spec.get("title", path.name)
            if not path.exists():
                lines.append(latex_escape(f"Missing figure: {path}"))
                continue
            rel_path = path.relative_to(output_path.parent) if path.is_relative_to(output_path.parent) else path
            lines.extend([
                r"\begin{figure}[H]",
                r"\centering",
                r"\includegraphics[width=0.92\linewidth]{\detokenize{" + str(rel_path) + "}}",
                r"\caption{" + latex_escape(title) + "}",
                r"\end{figure}",
            ])

    lines.append(r"\end{document}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=PROJECT_ROOT / "config" / "experiment_config.safe_parallel.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    config = load_json(args.config)
    report = config.get("report", {})
    output_path = args.output or resolve_path(report.get("tex_output", "output/reports/repro_compare.tex"), PROJECT_ROOT)
    if output_path is None:
        raise ValueError("Could not resolve report output path.")
    build_tex(config, output_path)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
