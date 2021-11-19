library(dplyr)
library(tidyr)
library(readr)

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

orpheum_cds <- read_csv(snakemake@input[['csv']][1]) %>%
  mutate(true_category = "Coding") %>%
  left_join(sam, orpheum_cds, by = "read_id") %>%
  select(read_id, jaccard_in_peptide_db, predicted_category = category, true_category, 
         predicted_translation_frame = translation_frame, true_translation_frame)

orpheum_noncds <- read_csv(snakemake@input[['csv']][2]) %>%
  mutate(true_category = "Non-coding") %>%
  mutate(true_translation_frame = NA) %>%
  select(read_id, jaccard_in_peptide_db, predicted_category = category, true_category, 
         predicted_translation_frame = translation_frame, true_translation_frame)

results <- bind_rows(orpheum_cds, orpheum_noncds)

results_filtered <- results %>%
  filter(predicted_category %in% c("Coding", "Non-coding")) %>%
  group_by(read_id) %>%
  arrange(desc(jaccard_in_peptide_db)) %>%
  slice(n = 1)

# table of FP/FN/TP/TN at different cutoffs -------------------------------

multi_class_rates_pct <- function(confusion_matrix) {
  true_positives  <- diag(confusion_matrix)
  false_positives <- colSums(confusion_matrix) - true_positives
  true_positives  <- c(true_positives, sum(true_positives))
  false_positives  <- c(false_positives, sum(false_positives))
  df <- data.frame(coding_true_positives  = (true_positives / sum(confusion_matrix)*100)[1], 
                   noncoding_true_positives = (true_positives / sum(confusion_matrix)*100)[2],
                   total_true_positives = (true_positives / sum(confusion_matrix)*100)[3],
                   coding_false_positives = (false_positives/sum(confusion_matrix)*100)[1], 
                   noncoding_false_positives = (false_positives/sum(confusion_matrix)*100)[2],
                   total_false_positives = (false_positives/sum(confusion_matrix)*100)[3],
                   row.names = 1)
  return(df)
}

accuracy_at_different_jaccard <- function(results, cutoff) {
  results <- results %>%
    mutate(new_predicted_category = ifelse(jaccard_in_peptide_db >= cutoff, "Coding", "Non-coding"))
  
  ct <- table(results$new_predicted_category, results$true_category)
  cm <- as.matrix(ct)
  df <- multi_class_rates_pct(cm)
  df$jaccard_cutoff <- cutoff
  return(df)
}

df_all_cutoffs <- data.frame()
for(cutoff in seq(from = .1, to = .6, by = .025)){
  df_cutoffs <- accuracy_at_different_jaccard(results_filtered, cutoff)
  df_all_cutoffs <- bind_rows(df_all_cutoffs, df_cutoffs)
}

df_all_cutoffs$accession <- snakemake@wildcards[["acc"]]

write_tsv(df_all_cutoffs, snakemake@output[['tsv']])
