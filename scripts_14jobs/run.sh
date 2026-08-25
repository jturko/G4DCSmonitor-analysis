#!/bin/bash
set -eu
MACRO_DIR="${1:?usage: ./run.sh macros/<dir>}"
TAG="$(basename "$MACRO_DIR")"
mkdir -p logs
sbatch --job-name="start_all_${TAG}" submit_start_all.sh "$MACRO_DIR"

