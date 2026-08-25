#!/bin/bash
set -eu
MACRO_DIR="${1:?usage: ./run.sh macros/<dir>}"
TAG="$(basename "$MACRO_DIR")"
mkdir -p logs

MANIFEST="$MACRO_DIR/manifest.txt"
# Put neutron jobs FIRST so the longest work starts immediately.
{ ls -1 "$MACRO_DIR"/run_baseFuel*neutron*.mac 2>/dev/null
  ls -1 "$MACRO_DIR"/run_baseFuel*.mac | grep -v neutron; } > "$MANIFEST"

NJOBS=$(wc -l < "$MANIFEST")
echo "Submitting $NJOBS jobs from $MANIFEST"
sbatch --job-name="g4_${TAG}" --array=0-$((NJOBS-1))%128 submit_macro.sh "$MANIFEST"

