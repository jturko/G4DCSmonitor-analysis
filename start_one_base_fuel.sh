#!/usr/bin/env bash
#
# Run the 5 source macros for ONE base fuel assembly. Each macro loads its
# ROOT file once and internally produces all 6 hexant outputs.
#
# Usage:   ./start_one_base_fuel.sh <base_fuel> <macro_dir>
# Example: ./start_one_base_fuel.sh 0 macros/auto
#
set -eu

if [ $# -ne 2 ]; then
    echo "Usage: $0 <base_fuel> <macro_dir>" >&2
    exit 1
fi

B="$1"
B_PAD="$(printf '%02d' "$B")"
MACRO_DIR="$2"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

shopt -s nullglob 2>/dev/null || true
MACS=( "$SCRIPT_DIR"/"$MACRO_DIR"/run_baseFuel${B_PAD}_*.mac )

if [ "${#MACS[@]}" -eq 0 ]; then
    echo "[err] no macros matching run_baseFuel${B_PAD}_*.mac in $SCRIPT_DIR/$MACRO_DIR" >&2
    echo "      (did you run ./generate_macros.sh first?)" >&2
    exit 4
fi

echo "[info] base fuel $B -> running ${#MACS[@]} source macros (each = 6 hexant runs):"
for MAC in "${MACS[@]}"; do echo "         $(basename "$MAC")"; done
echo "------------------------------------------------------------------"

for MAC in "${MACS[@]}"; do
    echo "[run] G4DCSmonitor $MAC"
    G4DCSmonitor "$MAC"
done

