library(rtracklayer)
library(plyranges)

gff <- rtracklayer::import(snakemake@input[['gff']])
pseudo <- gff %>% filter(pseudo == "true")

bed_noncds <- import(snakemake@input[['noncds_bed']])

noncds_pseudo <- subsetByOverlaps(bed_noncds, pseudo, type=c("any"))

noncds_coding_scores <- read_csv("outputs/orpheum/gtdb-rs202/protein-k10/GCA_000008025.1_noncds.coding_scores.csv") %>%
  filter(category == "Coding") %>%
  separate(read_id, into = c("tmp1", "tmp2", "seqnames", "ranges"), extra = "drop", remove = F, sep = ":") %>%
  select(-tmp1, -tmp2) %>%
  mutate(ranges = gsub("\\(.*", "", ranges)) %>%
  separate(ranges, into = c("start", "end"), sep = "-", remove = F) %>%
  mutate(start = as.numeric(start),
         end   = as.numeric(end))

noncds_predicted_as_coding_ranges<- as_granges(noncds_coding_scores)
noncds_predicted_as_coding_that_fall_in_pseudogenes <- subsetByOverlaps(noncds_predicted_as_coding_ranges, noncds_pseudo)

write_tsv(as.data.frame(noncds_predicted_as_coding_that_fall_in_pseudogenes), snakemake@output[['pseudo_tsv']])

summary_df <- data.frame(accession = snakemake@wildcards[['acc_w_gff']],
                         noncds_predicted_as_coding = length(unique(noncds_coding_scores$read_id)),
                         noncds_predicted_as_coding_in_pseudogenes = length(unique(noncds_predicted_as_coding_that_fall_in_pseudogenes$read_id)))
write_tsv(summary_df, snakemake@output[['pseudo_tab']])