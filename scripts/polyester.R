library(polyester)
library(Biostrings)
library(rtracklayer)
library(BSgenome)

# Use fasta of interest to simulate reads

fasta <- snakemake@input[['fasta']]
outdir <- snakemake@params[['outdir']]
num_reps <- 1
reads_per_seq <- 100

# build count matrix
fastaStrSet <- readDNAStringSet(fasta)
countmat = matrix(reads_per_seq, nrow=length(readDNAStringSet(fasta)), ncol=num_reps)
# simulate reads
simulate_experiment_countmat(fasta          = fasta, 
                             readmat        = countmat, 
                             reportCoverage = TRUE, 
                             readlen        = 150, 
                             paired         = FALSE, 
                             outdir         = outdir, 
                             gzip           = TRUE)
