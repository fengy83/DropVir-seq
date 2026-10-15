#!/usr/bin/env python3
"""Assign DropVir-seq reads to the predefined combinatorial barcode whitelist.

Methods statements implemented here
-----------------------------------
"We assign sequencing reads to virion-derived compartments using the bead-derived
barcode carried in Read 2, matching barcode sequences against the predefined
1,769,472-member whitelist and correcting sequences that differ from a unique
whitelist entry by one nucleotide. We discard reads without a valid barcode
assignment."

"Reads whose barcode differs from a unique whitelist entry by one nucleotide are
corrected; those equidistant from two entries fail the uniqueness test and are
discarded, so substitution errors cause loss rather than misassignment."

"The first 45 bases of Read 2 contain the bead-derived barcode construct and the
remaining 105 bases represent viral sequence."

"We first demultiplex reads according to the predefined i5/i7 index combinations
assigned to individual barcoding reactions, then group them by corrected droplet
barcode within each reaction, treating identical droplet-barcode sequences
occurring in different barcoding reactions as distinct identities."

Matching is therefore applied in two ordered passes, which is what the Methods
text describes:

  pass 1  exact whitelist match.  A read that matches the whitelist exactly is
          assigned outright; correction is not attempted.
  pass 2  one-nucleotide correction, applied only when pass 1 finds nothing, and
          only when the one-mismatch candidate is unique.

A read with more than one distinct candidate identity in the pass that resolves
it is discarded as ambiguous.

Barcode geometry
----------------
The whitelist is BC1 (3-6 nt) x BC2 (10 nt) x BC3 (9 nt) = 24 x 384 x 192 =
1,769,472 identities (Methods, "Fragmentation and end-repair ..."; table S3).
The Methods place the bead-derived construct in the first 45 bases of Read 2 but
do not give a 0-based start coordinate for BC1 inside that construct. The script
therefore searches the 45-nt construct for a contiguous BC1+BC2+BC3 identity
unless --barcode-start is supplied. Only an unambiguous whitelist identity is
accepted.

Read 2 is trimmed at a fixed offset of 45 bases, as stated in the Methods.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path


def op(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_fastq(path):
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip("\n")
            fh.readline()
            q = fh.readline().rstrip("\n")
            if not q:
                raise ValueError(f"Truncated FASTQ: {path}")
            yield h.rstrip("\n"), s, q


def read_list(path):
    vals = []
    with open(path) as fh:
        for line in fh:
            x = line.strip().upper()
            if x and not x.startswith("#"):
                vals.append(x)
    if not vals:
        raise ValueError(f"Empty barcode list: {path}")
    if len(vals) != len(set(vals)):
        raise ValueError(f"Duplicate sequences in {path}")
    return vals


def build_hamming1_index(values):
    """Exact set plus a wildcard index giving all entries within one substitution."""
    exact = set(values)
    patt = defaultdict(set)
    bylen = defaultdict(list)
    for v in values:
        bylen[len(v)].append(v)
        for i in range(len(v)):
            patt[(len(v), v[:i] + "*" + v[i + 1:])].add(v)
    return exact, patt, bylen


def candidates(obs, exact, patt):
    """Whitelist entries at Hamming distance 0 or 1 from obs, as (entry, distance)."""
    if obs in exact:
        return [(obs, 0)]
    hits = set()
    for i in range(len(obs)):
        hits.update(patt.get((len(obs), obs[:i] + "*" + obs[i + 1:]), ()))
    out = []
    for h in hits:
        d = sum(a != b for a, b in zip(obs, h))
        if d <= 1:
            out.append((h, d))
    return out


def segment_length_combinations(idxs):
    """Allowed (len BC1, len BC2, len BC3) triples.

    BC1 is 3-6 nt; BC2 and BC3 are fixed at 10 and 9 nt (Methods, table S3).
    """
    l1s = sorted(idxs[0][2])
    l2s = sorted(idxs[1][2])
    l3s = sorted(idxs[2][2])
    if len(l2s) != 1 or len(l3s) != 1:
        raise ValueError(
            "BC2 and BC3 must each have a single fixed length "
            f"(observed BC2 lengths {l2s}, BC3 lengths {l3s})"
        )
    return [(l1, l2s[0], l3s[0]) for l1 in l1s]


def collect_candidates(window, idxs, starts):
    """All whitelist identities within one substitution of the observed window.

    Returns {(bc1, bc2, bc3): (distance, start, end)} keeping the smallest
    distance when the same identity is reachable from more than one offset.
    """
    (e1, p1, _), (e2, p2, _), (e3, p3, _) = idxs
    combos = segment_length_combinations(idxs)
    found = {}
    for st in starts:
        for l1, l2, l3 in combos:
            end = st + l1 + l2 + l3
            if st < 0 or end > len(window):
                continue
            c1 = candidates(window[st:st + l1], e1, p1)
            if not c1:
                continue
            c2 = candidates(window[st + l1:st + l1 + l2], e2, p2)
            if not c2:
                continue
            c3 = candidates(window[st + l1 + l2:end], e3, p3)
            if not c3:
                continue
            for a, d1 in c1:
                for b, d2 in c2:
                    if d1 + d2 > 1:
                        continue
                    for c, d3 in c3:
                        d = d1 + d2 + d3
                        if d > 1:
                            continue
                        key = (a, b, c)
                        old = found.get(key)
                        if old is None or d < old[0]:
                            found[key] = (d, st, end)
    return found


def match_barcode(window, idxs, starts):
    """Exact match first; one-nucleotide correction only if no exact match exists.

    Returns (hit, status) with status in
    {"exact", "corrected", "ambiguous", "unassigned"}.
    """
    found = collect_candidates(window, idxs, starts)
    if not found:
        return None, "unassigned"

    exact = {k: v for k, v in found.items() if v[0] == 0}
    if exact:
        pool, status = exact, "exact"
    else:
        pool, status = found, "corrected"

    if len(pool) != 1:
        return None, "ambiguous"

    (a, b, c), (d, st, end) = next(iter(pool.items()))
    return (
        {"bc1": a, "bc2": b, "bc3": c, "barcode": a + b + c,
         "mismatches": d, "start": st, "end": end},
        status,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--bc1", required=True, help="24 first-segment sequences (table S3)")
    ap.add_argument("--bc2", required=True, help="384 second-segment sequences (table S3)")
    ap.add_argument("--bc3", required=True, help="192 third-segment sequences (table S3)")
    ap.add_argument("--reaction-id", required=True,
                    help="i5/i7-defined barcoding reaction identifier")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--barcode-window", type=int, default=45,
                    help="Read-2 barcode-construct length (Methods: 45 nt)")
    ap.add_argument("--barcode-start", type=int, default=None,
                    help="Known 0-based BC1 start inside the construct; "
                         "omit to search the window")
    ap.add_argument("--r2-viral-start", type=int, default=45,
                    help="First viral base of Read 2 (Methods: 45)")
    args = ap.parse_args()

    bc1, bc2, bc3 = (read_list(p) for p in (args.bc1, args.bc2, args.bc3))

    # The manuscript defines 24 x 384 x 192 = 1,769,472 identities, with
    # BC1 lengths of 3-6 nt, BC2 length 10 nt and BC3 length 9 nt. Fail
    # explicitly if a different design is supplied so the code cannot silently
    # drift from the reported barcode space.
    if len(bc1) != 24 or len(bc2) != 384 or len(bc3) != 192:
        raise ValueError(
            'Expected 24 BC1, 384 BC2 and 192 BC3 sequences '
            f'(observed {len(bc1)}, {len(bc2)}, {len(bc3)})'
        )
    if any(len(x) not in {3, 4, 5, 6} for x in bc1):
        raise ValueError('BC1 sequences must be 3-6 nt long')
    if any(len(x) != 10 for x in bc2):
        raise ValueError('BC2 sequences must be 10 nt long')
    if any(len(x) != 9 for x in bc3):
        raise ValueError('BC3 sequences must be 9 nt long')
    if len(bc1) * len(bc2) * len(bc3) != 1_769_472:
        raise ValueError('Barcode design must contain 1,769,472 identities')

    idxs = [build_hamming1_index(x) for x in (bc1, bc2, bc3)]
    lencombos = segment_length_combinations(idxs)
    shortest = min(l1 + l2 + l3 for l1, l2, l3 in lencombos)

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    stats = defaultdict(int)

    with gzip.open(out / "assigned_R1.fastq.gz", "wt") as r1out, \
         gzip.open(out / "assigned_R2.fastq.gz", "wt") as r2out, \
         gzip.open(out / "assignments.tsv.gz", "wt") as tout:
        tout.write("read_id\treaction_id\tbarcode\tbc1\tbc2\tbc3\t"
                   "mismatches\tbarcode_start\tbarcode_end\n")

        for n, (x1, x2) in enumerate(zip(read_fastq(args.r1), read_fastq(args.r2)), 1):
            h1, s1, q1 = x1
            h2, s2, q2 = x2
            id1 = h1.split()[0].lstrip("@")
            id2 = h2.split()[0].lstrip("@")
            if id1 != id2:
                raise ValueError(f"R1/R2 read-name mismatch at pair {n}: {id1} != {id2}")

            stats["total_pairs"] += 1
            win = s2[:args.barcode_window].upper()
            if args.barcode_start is not None:
                starts = [args.barcode_start]
            else:
                starts = range(max(0, len(win) - shortest) + 1)

            hit, status = match_barcode(win, idxs, starts)
            if hit is None:
                stats[status] += 1
                continue

            viral2 = s2[args.r2_viral_start:]
            viralq2 = q2[args.r2_viral_start:]
            if not viral2:
                stats["discarded_empty_r2_after_trim"] += 1
                continue

            stats[status] += 1
            stats["assigned"] += 1
            tag = (f" reaction={args.reaction_id} barcode={hit['barcode']}"
                   f" mismatches={hit['mismatches']}")
            r1out.write(f"@{id1}{tag}\n{s1}\n+\n{q1}\n")
            r2out.write(f"@{id2}{tag}\n{viral2}\n+\n{viralq2}\n")
            tout.write(f"{id1}\t{args.reaction_id}\t{hit['barcode']}\t{hit['bc1']}\t"
                       f"{hit['bc2']}\t{hit['bc3']}\t{hit['mismatches']}\t"
                       f"{hit['start']}\t{hit['end']}\n")

    with open(out / "barcode_assignment_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"reaction_id\t{args.reaction_id}\n")
        for k in ("total_pairs", "assigned", "exact", "corrected",
                  "ambiguous", "unassigned", "discarded_empty_r2_after_trim"):
            fh.write(f"{k}\t{stats[k]}\n")
    print(dict(stats), file=sys.stderr)


if __name__ == "__main__":
    main()
