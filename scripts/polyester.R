library(polyester)
library(Biostrings)
library(rtracklayer)
library(BSgenome)

# Use fasta of interest to simulate reads

fasta <- snakemake@input[['fasta']]
outdir <- snakemake@params[['outdir']]
num_reps <- as.integer(snakemake@params[["num_reps"]]) # 1
read_length <- as.integer(snakemake@params[["read_length"]]) # 150
paired <- snakemake@params[["simulate_paired"]] # true/false
reads_per_seq <- as.integer(snakemake@params[["num_reads_per_seq"]])

# build count matrix
fastaStrSet <- readDNAStringSet(fasta)
countmat = matrix(reads_per_seq, nrow=length(readDNAStringSet(fasta)), ncol=num_reps)
# simulate reads
simulate_experiment_countmat(fasta          = fasta, 
                             readmat        = countmat, 
                             reportCoverage = TRUE, 
                             readlen        = read_length, 
                             paired         = paired, 
                             outdir         = outdir, 
                             gzip           = TRUE)