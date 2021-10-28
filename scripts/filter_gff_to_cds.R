library(rtracklayer)

# import gff annotation files and filter to only coding domain sequences 
# (CDS), then rewrite in GFF format
gff_cds <- readGFF(snakemake@input[['gff']], filter = list(type=c("CDS")))
export(gff_cds, snakemake@output[['cds']])