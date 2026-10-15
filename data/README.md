# Barcode sequence inputs

The current preprocessing code uses the complete combinatorial barcode design described in the manuscript rather than the legacy subset barcode file.

Provide the barcode segments from table S3 as three plain-text files, one sequence per line:

- `BC1.txt`: 24 sequences, each 3-6 nt
- `BC2.txt`: 384 sequences, each 10 nt
- `BC3.txt`: 192 sequences, each 9 nt

Together these define 24 x 384 x 192 = 1,769,472 droplet-barcode identities.

`barcode3.txt` is retained only as a deprecation marker so that an older checkout does not silently keep using the previous subset whitelist. It must not be used as input to the revised pipeline.
