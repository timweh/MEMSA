# MEMSA - A MEM Extracting Multiple Sequence Aligner

MEMSA is a mutltiple sequence alignment (MSA) tool, which identifies maximum exact matches (MEMs) before applying a traditional MSA algorithm, in order to speed up the alignment process. It was developed to investigate the effects of this approximation on computation time and alignment accuracy and demonstrates that the preprocessing step can indeed positively impact the alignment of genomic sequences.

## Installation



## Usage

In order to use it, the MAFFT command line tools need to be preinstalled on the system, as MAFFT is called in order to perform alignment of the resulting subsequences.  

Also, slaMEM needs to be installed into the slaMEM directory.

The code can be executed by putting an arbitrary reference sequence in the reference.fa file and putting the other sequences to be aligned in the input.fa file. The generated alignment will be written into the alignment.fa file.
