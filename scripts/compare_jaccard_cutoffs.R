library(dplyr)
library(tidyr)
library(readr)

#sam <- read_tsv("outputs/cds_fasta_paladin/GCA_001884725.2_cds.sam", comment = "@", col_names = F)
sam <- read_tsv(snakemake@input[['sam']], comment = "@", col_names = F)
sam <- sam %>%
  select(X1) %>%
  mutate(read_id = gsub("^[0-9]+\\:[0-5]\\:[0-5]\\:", "", X1)) %>%
  separate(X1, into = c("seq_number", "orf_number", "relative_orf_number"), extra = "drop") %>%
  mutate(seq_number = as.numeric(seq_number) + 1) %>%
  mutate(orf_number_corrected = as.numeric(orf_number) + 1) %>%
  mutate(orf_number_corrected = ifelse(orf_number_corrected == 4, -1, orf_number_corrected)) %>%
  mutate(orf_number_corrected = ifelse(orf_number_corrected == 5, -2, orf_number_corrected)) %>%
  mutate(orf_number_corrected = ifelse(orf_number_corrected == 6, -3, orf_number_corrected)) %>%
  select(read_id, true_translation_frame = orf_number_corrected)

#orpheum_cds <- read_csv("outputs/orpheum/gtdb-rs202/protein-k10/GCA_001884725.2_cds.coding_scores.csv") %>%
orpheum_cds <- read_csv(snakemake@input[['csv']][1]) %>%
  left_join(sam, orpheum_cds, by = "read_id") %>%
  select(read_id, jaccard_in_peptide_db, predicted_category = category,
         predicted_translation_frame = translation_frame, true_translation_frame) %>%
  mutate(true_category = ifelse(true_translation_frame == predicted_translation_frame, "Coding", "Non-coding")) %>%
  select(read_id, jaccard_in_peptide_db, predicted_category, true_category, 
         predicted_translation_frame, true_translation_frame) 

orpheum_noncds <- read_csv("outputs/orpheum/gtdb-rs202/protein-k10/GCA_001884725.2_noncds.coding_scores.csv") %>%
  mutate(true_category = "Non-coding") %>%
  mutate(true_translation_frame = NA) %>%
  select(read_id, jaccard_in_peptide_db, predicted_category = category, true_category, 
         predicted_translation_frame = translation_frame, true_translation_frame)

results <- bind_rows(orpheum_cds, orpheum_noncds)

results_filtered <- results %>%
  filter(predicted_category %in% c("Coding", "Non-coding")) 


# table of FP/FN/TP/TN at different cutoffs -------------------------------

accuracy_at_different_jaccard <- function(results, cutoff) {
  results <- results %>%
    mutate(new_predicted_category = ifelse(jaccard_in_peptide_db >= cutoff, "Coding", "Non-coding"))
  
  ct <- table(results$true_category, results$new_predicted_category)
  cm <- as.matrix(ct)
  
  tp <- cm[1]
  tn <- cm[4]
  fp <- cm[2]
  fn <- cm[3]
  
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  
  youden_index = sensitivity + specificity - 1 
  
  df <- data.frame(tp, tn, fp, fn, sensitivity, specificity, youden_index)
  df$jaccard_cutoff <- cutoff
  return(df)
}

df_all_cutoffs <- data.frame()
for(cutoff in seq(from = 0.01, to = 1, by = .01)){
  df_cutoffs <- accuracy_at_different_jaccard(results_filtered, cutoff)
  df_all_cutoffs <- bind_rows(df_all_cutoffs, df_cutoffs)
}

df_all_cutoffs$accession <- snakemake@wildcards[["acc"]]

write_tsv(df_all_cutoffs, snakemake@output[['tsv']])