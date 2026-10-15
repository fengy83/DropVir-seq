#!/usr/bin/env python3
"""Reference assignment, purity and genome coverage in the defined community.

The filename is retained for compatibility with the original repository; the implementation now follows the manuscript reference-assignment, purity and coverage definitions.

Methods, "Reference assignment, purity and genome coverage in the defined
community":

  "We align defined-community assemblies to the T1, T4 and T7 reference genomes
   using minimap2 v2.28 with the -x asm5 preset and assign each vSAG to the
   reference receiving the greatest aligned length.  We calculate purity across
   all contigs associated with a barcode as the fraction of total alignable
   contig length assigned to the dominant phage reference, after merging
   overlapping query-aligned intervals, and classify a vSAG as pure only when
   100% of its alignable sequence is assigned to a single reference genome."

  "We calculate genome coverage as the union breadth of alignment to the
   assigned reference, that is, the ratio of reference bases covered by at least
   one sequence from that vSAG to the total number of bases in that reference
   genome, merging overlapping reference intervals before calculation."

  "We define high-confidence complete vSAGs as those that are pure and reach
   100% genome coverage."

Every alignment reported by minimap2 is retained; the defined-community Methods
do not specify suppression of secondary alignments for this analysis.
"""
import argparse
import csv
import os
import subprocess
from collections import defaultdict
from pathlib import Path


def fasta_lengths(path):
    d = {}
    name = None
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    d[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        d[name] = n
    return d


def merge(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = []
    s, e = intervals[0]
    for a, b in intervals[1:]:
        if a <= e:
            e = max(e, b)
        else:
            out.append((s, e))
            s, e = a, b
    out.append((s, e))
    return out


def ulen(intervals):
    return sum(b - a for a, b in merge(intervals))


def refmap(spec):
    """'label:seqid,label:seqid' -> {seqid: label}"""
    d = {}
    for x in spec.split(','):
        lab, seq = x.split(':', 1)
        d[seq] = lab
    return d


def paf_metrics(paf, rmap):
    """Per-reference merged query length, total alignable query length and
    per-reference target intervals."""
    q_by_ref = defaultdict(lambda: defaultdict(list))
    t_by_ref = defaultdict(list)
    all_q = defaultdict(list)
    with open(paf) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 12:
                continue
            q, qs, qe = p[0], int(p[2]), int(p[3])
            t, ts, te = p[5], int(p[7]), int(p[8])
            if t not in rmap:
                continue
            lab = rmap[t]
            q_by_ref[lab][q].append((qs, qe))
            all_q[q].append((qs, qe))
            t_by_ref[lab].append((ts, te))
    ref_bases = {lab: sum(ulen(v) for v in qd.values()) for lab, qd in q_by_ref.items()}
    total_alignable = sum(ulen(v) for v in all_q.values())
    return ref_bases, total_alignable, t_by_ref


def run_minimap2(contigs, refs, paf, threads):
    with open(paf, 'w') as out:
        subprocess.run(['minimap2', '-x', 'asm5', '-t', str(threads), refs, contigs],
                       stdout=out, check=True)


def main():
    a = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument('-i', '--input', required=True,
                   help='Directory tree of per-barcode assemblies')
    a.add_argument('-r', '--reference', required=True,
                   help='FASTA of the T1, T4 and T7 reference genomes')
    a.add_argument('-v', '--ref-ids', required=True,
                   help="'T1:NC_005833.1,T4:NC_000866.4,T7:NC_001604.1'")
    a.add_argument('-o', '--output', required=True)
    a.add_argument('-c', '--contig-filename', default='contigs.filtered.fasta')
    a.add_argument('-t', '--threads', type=int, default=4)
    args = a.parse_args()

    rmap = refmap(args.ref_ids)
    rlen = fasta_lengths(args.reference)
    labels = sorted(set(rmap.values()))
    out = Path(args.output)
    (out / 'paf').mkdir(parents=True, exist_ok=True)

    rows = []
    for root, _, files in os.walk(args.input):
        if args.contig_filename not in files:
            continue
        bid = os.path.basename(root)
        paf = str(out / 'paf' / f'{bid}.paf')
        run_minimap2(os.path.join(root, args.contig_filename),
                     args.reference, paf, args.threads)
        ref_bases, total, t_iv = paf_metrics(paf, rmap)

        if not ref_bases or total == 0:
            rows.append([bid, 'NA', 0, 0, '', 'Unassigned', ''] + [0] * len(labels))
            continue

        main_ref = max(ref_bases, key=ref_bases.get)
        dominant = ref_bases[main_ref]
        purity = dominant / total
        foreign = sum(v for k, v in ref_bases.items() if k != main_ref)
        pure = (foreign == 0)

        denom = sum(rlen.get(s, 0) for s, lab in rmap.items() if lab == main_ref)
        coverage = ulen(t_iv.get(main_ref, [])) / denom if denom else 0.0
        high_conf = pure and coverage >= 1.0 - 1e-12

        rows.append([
            bid, main_ref, total, dominant, f'{purity * 100:.6f}',
            'Pure' if pure else 'Impure', f'{coverage * 100:.6f}',
        ] + [ref_bases.get(x, 0) for x in labels] + [int(high_conf)])

    header = (['barcode', 'assigned_reference', 'total_alignable_bases',
               'dominant_aligned_bases', 'purity_percent', 'purity_class',
               'genome_coverage_percent']
              + [x + '_aligned_bases' for x in labels]
              + ['high_confidence_complete'])
    # Unassigned rows carry no coverage or completeness call; pad them to the
    # full width so no column is shifted.
    width = len(header)
    rows = [r + [''] * (width - len(r)) for r in rows]

    with open(out / 'vSAG_purity_coverage.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)
    print(f'{len(rows)} barcodes written to {out / "vSAG_purity_coverage.tsv"}')


if __name__ == '__main__':
    main()
