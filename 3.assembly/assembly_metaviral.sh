#!/usr/bin/env bash
# Methods: "We assemble reads assigned to each droplet barcode independently
# using SPAdes v4.2.0 in -metaviral mode with default k-mer selection, without
# pooling reads from different barcode-defined compartments."
#
# No -k is passed, so SPAdes applies its default k-mer selection.
set -euo pipefail

usage() { echo "Usage: $0 -i PER_BARCODE_DIR -o OUTDIR [-t THREADS_PER_JOB] [-j JOBS]" >&2; exit 2; }
THREADS=8; JOBS=4
while getopts 'i:o:t:j:h' x; do
  case "$x" in
    i) IN=$OPTARG ;; o) OUT=$OPTARG ;; t) THREADS=$OPTARG ;; j) JOBS=$OPTARG ;; *) usage ;;
  esac
done
[[ -n "${IN:-}" && -n "${OUT:-}" ]] || usage
mkdir -p "$OUT"

find "$IN" -type f \( -name '*_R1_001.fastq.gz' -o -name '*_R1_001.fastq' \) | sort > "$OUT/r1.list"

assemble_one() {
  r1=$1
  stem=$(basename "$r1"); stem=${stem%_R1_001.fastq.gz}; stem=${stem%_R1_001.fastq}
  r2=${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}
  [[ "$r2" == "$r1" ]] && r2=${r1/_R1_001.fastq/_R2_001.fastq}
  [[ -f "$r2" ]] || { echo "Missing mate: $r2" >&2; return 1; }
  od="$OUT/$stem"; mkdir -p "$od"
  spades.py -1 "$r1" -2 "$r2" --metaviral -o "$od" -t "$THREADS" > "$od/spades.log" 2>&1
}
export OUT THREADS
export -f assemble_one

xargs -a "$OUT/r1.list" -I{} -P "$JOBS" bash -c 'assemble_one "$1"' _ {}
echo "Per-barcode SPAdes v4.2.0 --metaviral assemblies written to $OUT"
