## Evaluating orpheum for open reading frame prediction from bacterial and archaeal short shotgun sequencing reads, using simulated reads

**Problem Statement**: Determine the correct open reading frame for bacterial and archaeal FASTQ sequencing reads to enable protein encodings for k-mers in the sequences.

**Approach**: Generate simulated reads from coding and non-coding portions of genomes and test whether orpheum can predict the correct open reading frame. Use genomes that are in the database (GTDB rs202) and not in the database (RefSeq, but added after the construction of GTDB rs202).

**Results**: In general, orpheum performed open reading frame prediction with high sensitivity and specificity. The results are analyzed in the `notebooks` directory, and written up in the repository github.com/dib-lab/2021-paper-metapangenomes.

## Getting started with this repository

The workflow is written in snakemake with software managed by conda. 
The workflow can be re-run with the following (snakemake command assumes a slurm cluster, and executes as a dry run first).
```
conda env create --name orph --file environment.yml
conda activate orph

snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=500000 --cluster "sbatch -t 10080 -J orph -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```
