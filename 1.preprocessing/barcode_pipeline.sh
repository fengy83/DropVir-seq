#!/usr/bin/env bash
# DropVir-seq preprocessing for one i5/i7-demultiplexed barcoding reaction.
#
# Materials and Methods:
#   1) quality-filter and adapter-trim paired-end reads with fastp;
#   2) assign reads using the predefined 1,769,472-member droplet-barcode
#      whitelist carried in Read 2, allowing correction only when a read differs
#      by one nucleotide from a unique whitelist entry;
#   3) discard reads without a valid barcode assignment;
#   4) group corrected droplet barcodes within each barcoding reaction.
#
# The Methods do not specify non-default fastp filtering/trimming parameters, so
# fastp is run with its defaults apart from I/O/report paths and thread count.
set -euo pipefail

if [[ $# -lt 8 ]]; then
  echo "Usage: $0 R1.fastq.gz R2.fastq.gz REACTION_ID BC1.txt BC2.txt BC3.txt OUTDIR THREADS [BARCODE_START]" >&2
  exit 2
fi

R1=$1
R2=$2
REACTION=$3
BC1=$4
BC2=$5
BC3=$6
OUT=$7
THREADS=$8
START=${9:-}
DIR="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUT/qc" "$OUT/assigned" "$OUT/per_barcode"

fastp \
  -i "$R1" -I "$R2" \
  -o "$OUT/qc/R1.clean.fastq.gz" -O "$OUT/qc/R2.clean.fastq.gz" \
  --json "$OUT/qc/fastp.json" --html "$OUT/qc/fastp.html" \
  -w "$THREADS"

CMD=(python "$DIR/extract_barcode.py"
  --r1 "$OUT/qc/R1.clean.fastq.gz"
  --r2 "$OUT/qc/R2.clean.fastq.gz"
  --bc1 "$BC1" --bc2 "$BC2" --bc3 "$BC3"
  --reaction-id "$REACTION"
  --outdir "$OUT/assigned"
  --barcode-window 45
  --r2-viral-start 45)

[[ -n "$START" ]] && CMD+=(--barcode-start "$START")
"${CMD[@]}"

python "$DIR/cluster_split.py" \
  --r1 "$OUT/assigned/assigned_R1.fastq.gz" \
  --r2 "$OUT/assigned/assigned_R2.fastq.gz" \
  --assignments "$OUT/assigned/assignments.tsv.gz" \
  --outdir "$OUT/per_barcode"

echo "Done: $OUT/per_barcode"
