#!/usr/bin/env python3
"""Length filter on per-barcode assemblies.

Methods: "We retain assembled contigs of at least 4,000 bp for downstream
analysis, applying the same length threshold to the defined community and to the
faecal catalogue."
"""
import argparse
import os


def records(path):
    h = None
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if h is not None:
                    yield h, ''.join(seq)
                h = line.rstrip('\n')
                seq = []
            else:
                seq.append(line.strip())
    if h is not None:
        yield h, ''.join(seq)


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('-i', '--input', required=True,
                   help='Directory tree of per-barcode assemblies')
    a.add_argument('-f', '--filename', default='contigs.fasta')
    a.add_argument('-o', '--output-name', default='contigs.filtered.fasta')
    a.add_argument('-l', '--length', type=int, default=4000, help='Methods: 4000 bp')
    args = a.parse_args()
    if args.output_name == args.filename:
        raise SystemExit('--output-name must differ from --filename')

    n_in = n_out = 0
    for root, _, files in os.walk(args.input):
        if args.filename not in files:
            continue
        with open(os.path.join(root, args.output_name), 'w') as out:
            for h, s in records(os.path.join(root, args.filename)):
                n_in += 1
                if len(s) >= args.length:
                    n_out += 1
                    out.write(h + '\n' + s + '\n')
    print(f'{n_out} of {n_in} contigs retained at >= {args.length} bp')


if __name__ == '__main__':
    main()
