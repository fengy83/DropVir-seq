# DropVir-seq

Analysis code for droplet-based particle-resolved viral genome sequencing.

This repository keeps the scope of the originally distributed DropVir-seq code. It does **not** attempt to add every analysis described in the manuscript. The code that is present has been aligned with the corresponding descriptions and parameters in the current Materials and Methods.

## Implemented workflow

```text
Paired-end reads from one predefined i5/i7 barcoding reaction
        |
        v
[1] Preprocessing
    fastp quality filtering / adapter trimming
    -> predefined combinatorial droplet-barcode assignment from Read 2
       (24 x 384 x 192 = 1,769,472 identities; unique 1-nt correction)
    -> reaction-aware per-barcode FASTQ
        |
        v
[2] Optional contig-to-reference alignment
    minimap2 v2.28, -x asm5
        |
        v
[3] Per-barcode assembly
    SPAdes v4.2.0 --metaviral, default k-mer selection
    -> retain contigs >= 4,000 bp
    -> CheckV v1.0.3 / checkv-db-v1.5
        |
        v
[4] Defined-community validation
    reference assignment / purity / genome coverage
    contig-level chimerism
```

## Directory structure

```text
DropVir-seq/
├── 1.preprocessing/
│   ├── barcode_pipeline.sh
│   ├── extract_barcode.py
│   └── cluster_split.py
├── 2.alignment/
│   └── align_minimap2.sh
├── 3.assembly/
│   ├── assembly_metaviral.sh
│   └── checkv_quality.sh
├── 4.analysis/
│   ├── quality_assessment/
│   │   ├── filter_contigs.py
│   │   └── map_and_count.py
│   └── mock_validation/
│       └── calculate_purity.py
├── data/
│   ├── README.md
│   └── barcode3.txt
├── environment.yml
├── requirements.txt
└── LICENSE
```

## Installation

```bash
conda env create -f environment.yml
conda activate dropvir-seq
```

The CheckV database used for the manuscript is `checkv-db-v1.5` and is supplied separately from the software environment.

## Barcode inputs

The manuscript defines 24 BC1 sequences (3-6 nt), 384 BC2 sequences (10 nt), and 192 BC3 sequences (9 nt), giving 1,769,472 combinatorial identities. The complete sequences are study data (table S3) and are not duplicated in this repository.

Prepare three plain-text files, one sequence per line:

```text
data/BC1.txt
data/BC2.txt
data/BC3.txt
```

The first 45 bases of Read 2 are treated as the bead-derived barcode construct; Read 2 sequence after base 45 is retained as viral sequence. If the exact 0-based BC1 start within the 45-nt construct is known for the processed library, it can be supplied as the optional final argument to `barcode_pipeline.sh`; otherwise the script searches the 45-nt construct and accepts only an unambiguous whitelist identity.

## Usage

### 1. Preprocessing

Run separately for each valid i5/i7 barcoding reaction. Identical corrected droplet-barcode sequences in different reactions remain distinct because the reaction identifier is retained.

```bash
bash 1.preprocessing/barcode_pipeline.sh \
  reaction01_R1.fastq.gz reaction01_R2.fastq.gz reaction_01 \
  data/BC1.txt data/BC2.txt data/BC3.txt \
  output/reaction_01 16
```

### 2. Per-barcode assembly

```bash
bash 3.assembly/assembly_metaviral.sh \
  -i output/reaction_01/per_barcode \
  -o assembly/reaction_01 \
  -t 8 -j 8
```

SPAdes is run in `--metaviral` mode without a user-specified `-k`, so the default k-mer selection is used.

### 3. Retain assembled contigs >=4,000 bp

```bash
python 4.analysis/quality_assessment/filter_contigs.py \
  -i assembly/reaction_01 \
  -f contigs.fasta \
  -o contigs.filtered.fasta \
  -l 4000
```

### 4. CheckV completeness

```bash
bash 3.assembly/checkv_quality.sh \
  assembly/reaction_01/combined_contigs.fasta \
  checkv_output \
  /path/to/checkv-db-v1.5 \
  16
```

### 5. Defined-community reference assignment, purity and coverage

```bash
python 4.analysis/mock_validation/calculate_purity.py \
  -i assembly/reaction_01 \
  -r references/T1_T4_T7.fasta \
  -v 'T1:NC_005833.1,T4:NC_000866.4,T7:NC_001604.1' \
  -o results/mock_purity
```

The script uses minimap2 v2.28 with `-x asm5`, assigns each vSAG to the reference with the greatest aligned length, calculates purity from merged query-aligned intervals across all contigs belonging to that barcode, and calculates genome coverage as union breadth on the assigned reference. A vSAG is classified as pure only when no alignable sequence is attributable to another phage reference; high-confidence complete vSAGs additionally require 100% reference coverage.

### 6. Defined-community chimerism

```bash
python 4.analysis/quality_assessment/map_and_count.py \
  -i assembly/reaction_01 \
  -r references/T1_T4_T7.fasta \
  -v 'T1:NC_005833.1,T4:NC_000866.4,T7:NC_001604.1' \
  -o results/chimeras.tsv
```

This uses BLAST+ v2.14.1 `blastn` with the default `megablast` task, E-value `1e-10`, and word size 28. A scored contig aligning to exactly two distinct phage references is called chimeric; no minimum per-parent alignment length or query-coverage threshold is imposed.

## Scope

The manuscript contains additional analyses that were not part of the originally distributed GitHub code. They are intentionally not added here. This repository is limited to the computational steps represented in the original codebase, with conflicting parameters or definitions corrected to match the current Materials and Methods.

## Authors

Hao Li, Yu Feng, Liuyang Zhao

## License

MIT License. See [LICENSE](LICENSE).
