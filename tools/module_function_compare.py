#!/usr/bin/env python3
from __future__ import annotations

import ast
import difflib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / 'src'
ORIG = ROOT / 'original codes and data'
OUT = ROOT / 'output' / 'module_function_compare.md'

DIRECT = {
    'src/iv/panel_iv.py': 'original codes and data/panel_iv.py',
    'src/iv_var/estimation.py': 'original codes and data/iv_var_estimation.py',
    'src/lmn_var/estimation.py': 'original codes and data/estimation.py',
    'src/utils/regression.py': 'original codes and data/regression.py',
    'src/model/solve.py': 'original codes and data/solve.py',
}


def text(p: Path) -> str:
    return p.read_text(encoding='utf-8', errors='ignore')


def defs(p: Path):
    src = text(p)
    try:
        t = ast.parse(src)
    except Exception:
        return []
    out = []
    for n in ast.walk(t):
        if isinstance(n, ast.ClassDef):
            out.append(('class', n.name, n.lineno))
        elif isinstance(n, ast.FunctionDef):
            out.append(('def', n.name, n.lineno))
    return sorted(out, key=lambda x: (x[0], x[1], x[2]))


def best_match(p: Path, cands):
    a = text(p)
    best = None
    for c in cands:
        s = difflib.SequenceMatcher(None, a, text(c)).ratio()
        if best is None or s > best[1]:
            best = (c, s)
    return best


def fmt_defs(ds):
    return ', '.join([f"{k} {n}@L{ln}" for k, n, ln in ds[:30]])


def main():
    src_files = sorted(SRC.rglob('*.py'))
    orig_py = sorted(ORIG.glob('*.py'))
    lines = ['# Module/Function Compare (Fresh)\n']

    lines.append('## A. Direct module mapping\n')
    for s_rel, o_rel in DIRECT.items():
        s = ROOT / s_rel
        o = ROOT / o_rel
        if not s.exists() or not o.exists():
            continue
        ss = text(s); oo = text(o)
        sim = difflib.SequenceMatcher(None, ss, oo).ratio()
        sdefs = defs(s); odefs = defs(o)
        sset = {(k,n) for k,n,_ in sdefs}
        oset = {(k,n) for k,n,_ in odefs}
        only_s = sorted(sset - oset)
        only_o = sorted(oset - sset)
        both = sorted(sset & oset)

        lines.append(f"### {s_rel} <-> {o_rel}")
        lines.append(f"- similarity: `{sim:.4f}`")
        lines.append(f"- lines(new/old): `{len(ss.splitlines())}/{len(oo.splitlines())}`")
        lines.append(f"- common symbols: `{len(both)}`")
        lines.append(f"- new-only symbols ({len(only_s)}): `" + ', '.join([f"{k} {n}" for k,n in only_s[:40]]) + '`')
        lines.append(f"- old-only symbols ({len(only_o)}): `" + ', '.join([f"{k} {n}" for k,n in only_o[:40]]) + '`')
        lines.append('')

    lines.append('## B. Every src file: best original match\n')
    lines.append('| src file | best old file | similarity | src symbols |')
    lines.append('|---|---|---:|---|')
    for s in src_files:
        b, sim = best_match(s, orig_py)
        sdefs = defs(s)
        lines.append(f"| {s.relative_to(ROOT)} | {b.relative_to(ROOT)} | {sim:.4f} | {len(sdefs)} |")

    lines.append('\n## C. Function inventory by module\n')
    for module in sorted([p for p in SRC.iterdir() if p.is_dir()]):
        lines.append(f"### module `{module.name}`")
        for f in sorted(module.glob('*.py')):
            if f.name == '__init__.py':
                continue
            d = defs(f)
            lines.append(f"- `{f.relative_to(ROOT)}` ({len(d)} symbols)")
            if d:
                lines.append(f"  - {fmt_defs(d)}")
        lines.append('')

    OUT.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f'Wrote {OUT}')


if __name__ == '__main__':
    main()
