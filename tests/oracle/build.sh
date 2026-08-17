#!/usr/bin/env bash
# Compiles the reference MD.cpp (golden oracle) exactly as the repo's own Makefile does.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT="${1:-$SCRIPT_DIR/MD_oracle.exe}"
g++ -O3 -I"$REPO_ROOT" -L"$REPO_ROOT" "$REPO_ROOT/MD.cpp" -o "$OUT"
echo "Built oracle binary: $OUT"
