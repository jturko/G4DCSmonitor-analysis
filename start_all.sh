#!/usr/bin/env bash
set -eu
MACRO_DIR="$1"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

BASE_FUELS=(0 1 2 3 5 6 7 8 13 14 15 22 23 32)
for B in "${BASE_FUELS[@]}"; do
    "$SCRIPT_DIR/start_one_base_fuel.sh" "$B" "$MACRO_DIR"
done

