#!/usr/bin/env bash
# Contig-to-reference alignment for the defined-community analysis.
# Materials and Methods: minimap2 v2.28 with the -x asm5 preset.
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 CONTIGS.fasta REFERENCE.fasta OUT.paf [THREADS]" >&2
  exit 2
fi

CONTIGS=$1
REF=$2
OUT=$3
THREADS=${4:-4}

minimap2 -x asm5 -t "$THREADS" "$REF" "$CONTIGS" > "$OUT"
echo "Wrote $OUT" >&2
