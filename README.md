# MEMSA - A MEM Extracting Multiple Sequence Aligner

MEMSA is a mutltiple sequence alignment (MSA) tool, which identifies maximum exact matches (MEMs) before applying a traditional MSA algorithm, in order to speed up the alignment process. It was developed to investigate the effects of this heuristic on computation time and alignment accuracy and demonstrates that the preprocessing step can indeed positively impact the alignment of genomic sequences.

## Requirements

This tool was developed for MacOS and Linux.
In order to built it, the gcc compiler needs to be installed.

## Manual

The tool can be executed by putting a single sequence in the reference FASTA-file and putting the all other sequences to be aligned in the input FASTA-file. The reference sequence can be picked arbitrarily from the dataset, as the choice of the reference does not affect the alignment. The generated alignment will be written into the output file.

### Install
```bash
./install.sh
```

The installation script downloads the required dependencies [slaMEM](https://github.com/fjdf/slaMEM) and [MAFFT](https://mafft.cbrc.jp/alignment/software/) and builds an executable from the source code.

### Usage
```bash
./memsa <options>
```

To run MEMSA for the provided example files and default parameters, just run `./memsa`

#### Options:
- `-s`   : minimum seed length (default=20)
- `-g`   : maximum merge gap (default=1)
- `-r`   : reference file name (default="reference.fa")
- `-i`   : input file name (default="input.fa")
- `-o`   : output file name (default="alignment.fa")

#### Example:
```bash
./memsa
./memsa -s 50 -g 0
./memsa -r ref.fa -i sequences.fa -o result.fa
```