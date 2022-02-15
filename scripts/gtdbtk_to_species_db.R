library(dplyr)
library(readr)
library(tidyr)

gtdbtk <- read_tsv(snakemake@input[['gtdbtk_arc']], col_types = "cccdcddcddddccccdc", na = "N/A")
gtdbtk <- read_tsv(snakemake@input[['gtdbtk_bac']], col_types = "cccdcddcddddccccdc", na = "N/A") %>%
  bind_rows(gtdbtk) %>%
  mutate(user_genome = gsub("_genomic", "", user_genome)) %>%
  mutate(accession_minus_prefix = gsub("GCA_", "", user_genome)) %>%
  mutate(accession_minus_prefix = gsub("\\.[1-9].*", "", accession_minus_prefix)) %>%
  separate(classification, into=c("superkingdom", "phylum", "order", "class", "family", "genus", "species"), sep = ";")

gtdbtk_assigned_species <- gtdbtk %>%
  filter(species != "s__")

nrow(gtdbtk_assigned_species)

gtdbtk_assigned_species <- gtdbtk_assigned_species %>%
  select(accession = user_genome, species) %>%
  mutate(species = gsub(" ", "_", species))

write.csv(gtdbtk_assigned_species, snakemake@output[['csv']])
# use this to link in dbs to inputs/orpheum_index
# cat(paste0("/home/ntpierce/2021-rank-compare/output.rank-compare/sourmash-nodegraph/species/gtdb-rs202.", gtdbtk_assigned_species$species, ".protein-k10.nodegraph"))
# note that this db is missing for now: gtdb-rs202.s__Natrinema_sp013456555.protein-k10.nodegraph
