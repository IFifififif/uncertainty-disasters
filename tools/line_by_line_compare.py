#!/usr/bin/env python3
from __future__ import annotations

import ast
import difflib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
ORIG = ROOT / "original codes and data"
OUT_MD = ROOT / "output" / "code_line_by_line_audit.md"
OUT_JSON = ROOT / "output" / "code_line_by_line_audit.json"

DIRECT_MAP = {
    "src/iv/panel_iv.py": "original codes and data/panel_iv.py",
    "src/iv_var/estimation.py": "original codes and data/iv_var_estimation.py",
    "src/lmn_var/estimation.py": "original codes and data/estimation.py",
    "src/utils/regression.py": "original codes and data/regression.py",
    "run_all.py": "original codes and data/run_all.py",
    "src/model/solve.py": "original codes and data/solve.py",
}


@dataclass
class FileCmp:
    new_file: str
    old_file: str
    new_lines: int
    old_lines: int
    similarity: float
    new_only_defs: List[str]
    old_only_defs: List[str]
    common_defs: int
    diff_preview: List[str]


def read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="ignore")


def list_defs(py_text: str) -> List[str]:
    try:
        tree = ast.parse(py_text)
    except Exception:
        return []
    names = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.append(node.name)
    return sorted(set(names))


def unified_preview(a_text: str, b_text: str, a_name: str, b_name: str, max_lines: int = 220) -> List[str]:
    ud = list(
        difflib.unified_diff(
            a_text.splitlines(),
            b_text.splitlines(),
            fromfile=a_name,
            tofile=b_name,
            lineterm="",
            n=2,
        )
    )
    if len(ud) > max_lines:
        return ud[:max_lines] + [f"... ({len(ud)-max_lines} more diff lines truncated)"]
    return ud


def compare_pair(new_file: Path, old_file: Path) -> FileCmp:
    new_text = read_text(new_file)
    old_text = read_text(old_file)
    sim = difflib.SequenceMatcher(None, new_text, old_text).ratio()

    new_defs = set(list_defs(new_text))
    old_defs = set(list_defs(old_text))

    return FileCmp(
        new_file=str(new_file.relative_to(ROOT)),
        old_file=str(old_file.relative_to(ROOT)),
        new_lines=len(new_text.splitlines()),
        old_lines=len(old_text.splitlines()),
        similarity=sim,
        new_only_defs=sorted(new_defs - old_defs),
        old_only_defs=sorted(old_defs - new_defs),
        common_defs=len(new_defs & old_defs),
        diff_preview=unified_preview(old_text, new_text, str(old_file.relative_to(ROOT)), str(new_file.relative_to(ROOT))),
    )


def best_match(new_file: Path, candidates: List[Path]) -> Dict:
    txt = read_text(new_file)
    best = None
    for cand in candidates:
        ctxt = read_text(cand)
        sim = difflib.SequenceMatcher(None, txt, ctxt).ratio()
        if best is None or sim > best["similarity"]:
            best = {"old_file": str(cand.relative_to(ROOT)), "similarity": sim}
    return best


def main() -> None:
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)

    src_py = sorted([p for p in SRC.rglob("*.py")])
    orig_py = sorted([p for p in ORIG.glob("*.py")])

    # direct mapped detailed comparisons
    detailed: List[FileCmp] = []
    for n_rel, o_rel in DIRECT_MAP.items():
        n_path = ROOT / n_rel
        o_path = ROOT / o_rel
        if n_path.exists() and o_path.exists():
            detailed.append(compare_pair(n_path, o_path))

    # whole repo best-match view
    best = []
    for nf in src_py:
        b = best_match(nf, orig_py)
        best.append({
            "new_file": str(nf.relative_to(ROOT)),
            "best_old_file": b["old_file"],
            "similarity": b["similarity"],
        })

    best = sorted(best, key=lambda x: x["new_file"])

    payload = {
        "direct_map": [c.__dict__ for c in detailed],
        "best_match_all_src_py": best,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("# Code Line-by-Line Audit (Fresh)\n")
    lines.append("## 1) Direct Mapped Files (Detailed)\n")
    for c in detailed:
        lines.append(f"### {c.new_file}  vs  {c.old_file}\n")
        lines.append(f"- Similarity: `{c.similarity:.4f}`")
        lines.append(f"- Line count (new/old): `{c.new_lines}/{c.old_lines}`")
        lines.append(f"- Common defs: `{c.common_defs}`")
        lines.append(f"- New-only defs ({len(c.new_only_defs)}): `{', '.join(c.new_only_defs[:20])}`")
        lines.append(f"- Old-only defs ({len(c.old_only_defs)}): `{', '.join(c.old_only_defs[:20])}`")
        lines.append("- Diff preview:")
        lines.append("```diff")
        lines.extend(c.diff_preview)
        lines.append("```\n")

    lines.append("## 2) Best Match for Every src/*.py\n")
    lines.append("| new_file | best_old_file | similarity |")
    lines.append("|---|---|---:|")
    for row in best:
        lines.append(f"| {row['new_file']} | {row['best_old_file']} | {row['similarity']:.4f} |")

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_MD}")
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
