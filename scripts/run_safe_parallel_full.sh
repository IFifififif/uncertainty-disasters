#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

CONFIG="config/experiment_config.safe_parallel.json"
WITH_REPORT=0
DRY_RUN=0
SELECTED_GROUPS=""

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_safe_parallel_full.sh [options]

Options:
  --config PATH       JSON config to use. Default: config/experiment_config.safe_parallel.json
  --groups A,B        Comma-separated run_groups to execute. Default: all enabled groups in config.
  --with-report       Generate the configured LaTeX report after runs finish.
  --dry-run           Print resolved groups and commands without executing run_all.py.
  -h, --help          Show this help.

The config controls:
  - defaults.*.data_path for user datasets
  - jobs[].modules and job-level parameter overrides
  - run_groups[] for parallel/sequential execution phases
  - report.* for generated .tex content, labels, tables, figures, comparisons
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --config)
      CONFIG="$2"
      shift 2
      ;;
    --groups)
      SELECTED_GROUPS="$2"
      shift 2
      ;;
    --with-report)
      WITH_REPORT=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [ ! -d ".venv" ]; then
  echo "Error: .venv not found at $ROOT_DIR/.venv" >&2
  echo "Create it with: python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt" >&2
  exit 1
fi

if [ ! -f "$CONFIG" ]; then
  echo "Error: config not found: $CONFIG" >&2
  exit 1
fi

source .venv/bin/activate
mkdir -p output/test_logs output/tmp output/reports

PLAN_JSON="output/tmp/safe_parallel_plan.json"
AGG_SUMMARY="output/tmp/run_summary_aggregate.json"
printf '[]\n' > "$AGG_SUMMARY"
python - "$CONFIG" "$SELECTED_GROUPS" "$PLAN_JSON" <<'PY'
import json
import sys
from pathlib import Path

config_path = Path(sys.argv[1])
selected = {x.strip() for x in sys.argv[2].split(",") if x.strip()}
plan_path = Path(sys.argv[3])

cfg = json.loads(config_path.read_text(encoding="utf-8"))
groups = cfg.get("run_groups", [])
if not groups:
    raise SystemExit("Config must define at least one run_groups entry.")

chosen = []
for group in groups:
    name = group.get("name")
    if selected:
        if name in selected:
            chosen.append(group)
    elif group.get("enabled", True):
        chosen.append(group)

if not chosen:
    raise SystemExit("No run groups selected. Check --groups or run_groups[].enabled.")

known_jobs = {job.get("name") for job in cfg.get("jobs", [])}
for group in chosen:
    missing = [job for job in group.get("jobs", []) if job not in known_jobs]
    if missing:
        raise SystemExit(f"Run group {group.get('name')} references unknown jobs: {missing}")

plan_path.parent.mkdir(parents=True, exist_ok=True)
plan_path.write_text(json.dumps(chosen, ensure_ascii=False, indent=2), encoding="utf-8")
print(plan_path)
PY

python - "$PLAN_JSON" <<'PY'
import json
import sys
from pathlib import Path

groups = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
print("Selected run groups:")
for group in groups:
    mode = "sequential" if group.get("sequential", False) else "parallel"
    print(f"  - {group.get('name')} ({mode}): {', '.join(group.get('jobs', []))}")
PY

idx=0
while IFS= read -r group_name; do
  idx=$((idx + 1))
  GROUP_CONFIG="output/tmp/config_${group_name}.json"
  LOG_PATH="output/test_logs/${idx}_${group_name}.log"

  python - "$CONFIG" "$PLAN_JSON" "$group_name" "$GROUP_CONFIG" <<'PY'
import json
import sys
from pathlib import Path

config_path = Path(sys.argv[1])
plan_path = Path(sys.argv[2])
group_name = sys.argv[3]
out_path = Path(sys.argv[4])

cfg = json.loads(config_path.read_text(encoding="utf-8"))
groups = json.loads(plan_path.read_text(encoding="utf-8"))
group = next(g for g in groups if g.get("name") == group_name)
enabled_jobs = set(group.get("jobs", []))

for job in cfg.get("jobs", []):
    job["enabled"] = job.get("name") in enabled_jobs

cfg["parallel"] = dict(cfg.get("parallel", {}))
if group.get("sequential", False):
    cfg["parallel"]["enabled"] = False

out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")
print(out_path)
PY

  SEQ_FLAG=""
  if python - "$PLAN_JSON" "$group_name" <<'PY'
import json
import sys
from pathlib import Path
groups = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
group = next(g for g in groups if g.get("name") == sys.argv[2])
raise SystemExit(0 if group.get("sequential", False) else 1)
PY
  then
    SEQ_FLAG="--sequential"
  fi

  echo "[$idx] Running group: $group_name"
  echo "    config: $GROUP_CONFIG"
  echo "    log:    $LOG_PATH"
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "    dry-run command: python run_all.py --config $GROUP_CONFIG $SEQ_FLAG"
  else
    python run_all.py --config "$GROUP_CONFIG" $SEQ_FLAG 2>&1 | tee "$LOG_PATH"
    python - "$AGG_SUMMARY" "output/runs/run_summary.json" <<'PY'
import json
import sys
from pathlib import Path

agg_path = Path(sys.argv[1])
summary_path = Path(sys.argv[2])
agg = json.loads(agg_path.read_text(encoding="utf-8"))
latest = json.loads(summary_path.read_text(encoding="utf-8"))
agg.extend(latest)
agg_path.write_text(json.dumps(agg, ensure_ascii=False, indent=2), encoding="utf-8")
summary_path.write_text(json.dumps(agg, ensure_ascii=False, indent=2), encoding="utf-8")
PY
  fi
done < <(python - "$PLAN_JSON" <<'PY'
import json
import sys
from pathlib import Path
for group in json.loads(Path(sys.argv[1]).read_text(encoding="utf-8")):
    print(group.get("name"))
PY
)

if [ "$WITH_REPORT" -eq 1 ]; then
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "dry-run command: python tools/build_repro_compare_tex.py --config $CONFIG"
  else
    python tools/build_repro_compare_tex.py --config "$CONFIG"
  fi
fi

echo "Done."
