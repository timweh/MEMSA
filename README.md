# MEMSA - A MEM Extracting Multiple Sequence Aligner

MEMSA is a mutltiple sequence alignment (MSA) tool, which identifies maximum exact matches (MEMs) before applying a traditional MSA algorithm, in order to speed up the alignment process. It was developed to investigate the effects of this heuristic on computation time and alignment accuracy and demonstrates that the preprocessing step can indeed positively impact the alignment of genomic sequences.

## Requirements

This tool was developed for MacOS and Linux.
In order to built it, the gcc compiler needs to be installed.

## Installation

In order for the tool to be used, [slaMEM](https://github.com/fjdf/slaMEM) needs to be installed into the slaMEM directory and [MAFFT](https://mafft.cbrc.jp/alignment/software/) into the mafft directory. In order to install independencies and create all required subdirectories, the setup.sh script needs to be called before running MEMSA for the first time.

## Usage

The tool takes two fasta files as input. One is the reference file and must contain a single sequence. All other sequences to be aligned go into the input file. The reference sequence can be picked arbitrarily from the dataset, as the choice of the reference does not affect the alignment. 

### Input

### Output

### Command Line Arguments

The code can be executed by putting an arbitrary reference sequence in the reference.fa file and putting the other sequences to be aligned in the input.fa file. The generated alignment will be written into the alignment.fa file.

To run MEMSA for the provided example files and default parameters, just run `./memsa`

### Manual

#### Install
```bash
./install.sh
```

#### Usage
```bash
./memsa <options>
```

##### Options:
- `s`   : minimum seed length (default=20)
- `g`   : maximum merge gap (default=1)
- `r`   : reference file name (default="reference.fa")
- `i`   : input file name (default="input.fa")
- `o`   : output file name (default="alignment.fa")

##### Example:
```bash
./memsa
./memsa -s 50 -g 0
./memsa -r ref.fa -i sequences.fa -o result.fa
```