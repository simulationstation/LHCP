#!/usr/bin/env bash
set -euo pipefail

wdalpha holistic-run \
  --target G191-B2B \
  --data-root ./LHCP/ \
  --out ./LHCP/results/holistic_smoke/ \
  --atomic ./LHCP/data/atomic/analysis_lines_clean.csv \
  --config scripts/holistic_smoke.yaml
