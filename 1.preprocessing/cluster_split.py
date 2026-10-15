#!/usr/bin/env python3
"""Split assigned reads into reaction-aware per-barcode FASTQ files.

The filename is retained for compatibility with the original repository. No
sequence clustering is performed: the manuscript uses corrected predefined
whitelist identities rather than Starcode clusters.

Methods: "... then group them by corrected droplet barcode within each reaction,
treating identical droplet-barcode sequences occurring in different barcoding
reactions as distinct identities."

Output files are therefore keyed by reaction_id and barcode jointly.  No minimum
read count is applied, because the Methods define no such filter.
"""
import argparse
import gzip
import re
from collections import Counter, OrderedDict
from pathlib import Path


def rfq(path):
    with gzip.open(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                return
            yield h, f.readline(), f.readline(), f.readline()


def read_assign(path):
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {x: i for i, x in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            yield p[ix["read_id"]], p[ix["reaction_id"]], p[ix["barcode"]]


def clean(x):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", x)


class HandleCache:
    """Bounded pool of gzip handles.

    A key is truncated the first time it is written and appended to thereafter,
    so re-running the script does not concatenate onto a previous run's output.
    """

    def __init__(self, limit):
        self.limit = limit
        self.open = OrderedDict()
        self.started = set()

    def get(self, key, path):
        if key in self.open:
            h = self.open.pop(key)
            self.open[key] = h
            return h
        if len(self.open) >= self.limit:
            _, old = self.open.popitem(last=False)
            old.close()
        mode = "at" if key in self.started else "wt"
        self.started.add(key)
        h = gzip.open(path, mode)
        self.open[key] = h
        return h

    def close(self):
        for h in self.open.values():
            h.close()
        self.open.clear()


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--r1", required=True)
    a.add_argument("--r2", required=True)
    a.add_argument("--assignments", required=True)
    a.add_argument("--outdir", required=True)
    a.add_argument("--max-open", type=int, default=128,
                   help="Maximum simultaneously open output handles")
    args = a.parse_args()

    counts = Counter((r, b) for _, r, b in read_assign(args.assignments))
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    h1 = HandleCache(args.max_open)
    h2 = HandleCache(args.max_open)

    for fq1, fq2, rec in zip(rfq(args.r1), rfq(args.r2), read_assign(args.assignments)):
        rid, reaction, bc = rec
        id1 = fq1[0].split()[0].lstrip("@")
        id2 = fq2[0].split()[0].lstrip("@")
        if id1 != rid or id2 != rid:
            raise ValueError("FASTQ/assignment order mismatch")
        stem = f"{clean(reaction)}__{bc}_S1_L001"
        h1.get(stem, str(out / f"{stem}_R1_001.fastq.gz")).writelines(fq1)
        h2.get(stem, str(out / f"{stem}_R2_001.fastq.gz")).writelines(fq2)

    h1.close()
    h2.close()

    with open(out / "barcode_counts.tsv", "w") as f:
        f.write("reaction_id\tbarcode\tread_pairs\n")
        for (r, b), n in sorted(counts.items()):
            f.write(f"{r}\t{b}\t{n}\n")


if __name__ == "__main__":
    main()
