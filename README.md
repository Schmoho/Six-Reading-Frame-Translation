# Six Reading Frame Translation Script

## Usage Notes

The script can be called with `python translate.py -code path/to/aa-code -intervals path/to/intervals -sequences path/to/fasta`.
Assume relative paths to the point of invocation.

The implementation does not make any effort to ensure viable input or fail gracefully.

No external packages are required, however `conda` environment and `pip` requirements.txt file are provided to suggest a suitable python version.

## What this does:

Extracts and translates genes with defined coordinates from their chromosome sequences.

Three test files are provided:

* **`sequences.fasta`** A FASTA file with two chromosome sequences.
* **`intervals.gff`** A GFF file with a set of mRNA coordinates for the sequences in `sequences.fasta`.
* **`standard_code.txt`** A tab delimited table mapping codons to single letter amino acids.

1. Extracts the mRNA sequences from `sequences.fasta` given the coordinates in `intervals.gff`,
2. Translates all six reading frames of the mRNA sequence given the genetic code in `standard_code.txt`,
