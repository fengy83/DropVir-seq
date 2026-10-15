#!/usr/bin/env bash
# Methods: "We estimate genome completeness independently using CheckV v1.0.3
# with checkv-db-v1.5."
#
# Downstream, the faecal catalogue retains vSAGs with estimated completeness of
# at least 90% as near-complete.  Apply that filter to the CheckV
# quality_summary.tsv produced here; no additional CheckV options are specified
# in the Methods, so none are passed.
set -euo pipefail
[[ $# -ge 3 ]] || { echo "Usage: $0 INPUT_FASTA OUTDIR CHECKV_DB [THREADS]" >&2; exit 2; }
IN=$1; OUT=$2; DB=$3; T=${4:-16}
checkv end_to_end "$IN" "$OUT" -d "$DB" -t "$T"
