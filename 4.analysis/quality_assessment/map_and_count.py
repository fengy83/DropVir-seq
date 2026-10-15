#!/usr/bin/env python3
"""Contig-level chimerism in the defined community.

The filename is retained for compatibility with the original repository; the implementation now follows the manuscript chimerism definition.

Methods, "Sequence fidelity, amplification uniformity and chimerism":

  "We align each scored contig against the concatenated T1, T4 and T7 references
   using blastn from BLAST+ v2.14.1 with an expectation-value threshold of
   1 x 10-10, the default megablast task and a word size of 28, and classify a
   contig producing alignments to exactly two distinct phage references as
   chimeric, imposing no minimum per-parent alignment length or query-coverage
   threshold.  We evaluate chimera frequency per contig."

The Methods define the chimeric class as "exactly two distinct phage
references".  Contigs aligning to three references therefore fall outside the
definition; they are reported in their own class and counted separately in the
summary rather than silently merged into the non-chimeric count.

A scored contig is a contig producing at least one alignment to a reference.
The summary also reports the number of contigs examined, so both possible
denominators are visible.
"""
import argparse
import csv
import os
import subprocess
from collections import defaultdict


def parse_map(spec):
    d = {}
    for x in spec.split(','):
        lab, seq = x.split(':', 1)
        d[seq] = lab
    return d


def count_contigs(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                n += 1
    return n


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
    a.add_argument('-f', '--filename', default='contigs.filtered.fasta')
    a.add_argument('-t', '--threads', type=int, default=1)
    args = a.parse_args()

    rm = parse_map(args.ref_ids)
    rows = []
    n_examined = 0

    for root, _, files in os.walk(args.input):
        if args.filename not in files:
            continue
        bid = os.path.basename(root)
        query = os.path.join(root, args.filename)
        n_examined += count_contigs(query)
        cmd = ['blastn', '-task', 'megablast', '-query', query,
               '-subject', args.reference, '-evalue', '1e-10',
               '-word_size', '28', '-outfmt', '6 qseqid sseqid',
               '-num_threads', str(args.threads)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        hits = defaultdict(set)
        for line in res.splitlines():
            q, s = line.split('\t')[:2]
            if s in rm:
                hits[q].add(rm[s])
        for q, hs in hits.items():
            n = len(hs)
            if n == 2:
                call = 'Chimeric'
            elif n == 1:
                call = 'Non_chimeric'
            else:
                call = 'Multi_reference'
            rows.append([bid, q, ','.join(sorted(hs)), n, call])

    with open(args.output, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['barcode', 'contig', 'reference_parents',
                    'n_distinct_references', 'chimera_call'])
        w.writerows(rows)

    scored = len(rows)
    chim = sum(r[-1] == 'Chimeric' for r in rows)
    multi = sum(r[-1] == 'Multi_reference' for r in rows)
    with open(args.output + '.summary.txt', 'w') as f:
        f.write(f'contigs_examined\t{n_examined}\n')
        f.write(f'scored_contigs\t{scored}\n')
        f.write(f'chimeric_contigs\t{chim}\n')
        f.write(f'contigs_with_more_than_two_references\t{multi}\n')
        f.write(f'frequency_over_scored_contigs\t{chim / scored if scored else 0:.8f}\n')


if __name__ == '__main__':
    main()
