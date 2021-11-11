library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
setwd('~/github/2021-orpheum-sim/')

sam <- read_tsv("sandbox/test_paladin_for_ORF_detection/tmp.sam", comment = "@", col_names = F)
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

orpheum_cds <- read_csv("outputs/orpheum/gtdb-rs202/protein-k10/GCA_000008025.1_cds.coding_scores.csv") %>%
  mutate(true_category = "Coding") %>%
  left_join(sam, orpheum_cds, by = "read_id") %>%
  select(read_id, jaccard_in_peptide_db, predicted_category = category, true_category, 
         predicted_translation_frame = translation_frame, true_translation_frame)

orpheum_noncds <- read_csv("outputs/orpheum/gtdb-rs202/protein-k10/GCA_000008025.1_noncds.coding_scores.csv") %>%
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




# plot results ------------------------------------------------------------
ggplot(results_filtered, aes(x = jaccard_in_peptide_db, fill = true_category)) +
  geom_histogram(bins = 100) +
  theme_minimal() +
  scale_y_sqrt()


# table of FP/FN/TP/TN at different cutoffs -------------------------------

multi_class_rates <- function(confusion_matrix) {
  true_positives  <- diag(confusion_matrix)
  false_positives <- colSums(confusion_matrix) - true_positives
  false_negatives <- rowSums(confusion_matrix) - true_positives
  true_negatives  <- sum(confusion_matrix) - true_positives -
    false_positives - false_negatives
  return(data.frame(true_positives, false_positives, true_negatives,
                    false_negatives, row.names = names(true_positives)))
}

multi_class_rates_pct <- function(confusion_matrix) {
  true_positives  <- diag(confusion_matrix)
  false_positives <- colSums(confusion_matrix) - true_positives
  true_positives  <- c(true_positives, sum(true_positives))
  false_positives  <- c(false_positives, sum(false_positives))
  return(data.frame(true_positives  = true_positives /sum(confusion_matrix)*100, 
                    false_positives = false_positives/sum(confusion_matrix)*100, 
                    row.names = c("Coding", "Non-coding", "Total")))
}

accuracy_at_different_jaccard <- function(results, cutoff) {
  results <- results %>%
    mutate(new_predicted_category = ifelse(jaccard_in_peptide_db >= cutoff, "Coding", "Non-coding"))
  
  ct <- table(results$new_predicted_category, results$true_category)
  cm <- as.matrix(ct)
  multi_class_rates_pct(cm)
}

accuracy_at_different_jaccard(results_filtered, .5)
accuracy_at_different_jaccard(results_filtered, .4)
accuracy_at_different_jaccard(results_filtered, .35)
accuracy_at_different_jaccard(results_filtered, .3)
accuracy_at_different_jaccard(results_filtered, .2)
accuracy_at_different_jaccard(results_filtered, .1)

# old plots ---------------------------------------------------------------


# only plot outliers as points in the plot
# results_filtered <- results_filtered %>%
#   group_by(category) %>%
#   mutate(outlier = jaccard_in_peptide_db > quantile(jaccard_in_peptide_db, .75) + 1.50*IQR(jaccard_in_peptide_db)) %>%
#   mutate(outlier = ifelse(outlier == TRUE, TRUE, jaccard_in_peptide_db < quantile(jaccard_in_peptide_db, .25) - 1.50*IQR(jaccard_in_peptide_db))) %>%
#   ungroup()
# Boxplot of jaccard similarity scores for all six frame translations of each read 
# when compared against a database of all 10-mers in GTDB rs202. Scores are separated 
# by accuracy of the translation frame of the read. Whiskers indicate the largest 
# and smallest value within 1.5 times the interquartile range below the 25th percentile.
# Outliers are jittered. No matter what cutoff is used, some false positives and
# false negatives will occur. We suggest using 0.4 as a jaccard similarity cutoff.

# ggplot(results_filtered, aes(x = category, y = jaccard_in_peptide_db)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(data = function(x) dplyr::filter_(x, ~ outlier), position = 'jitter',
#              alpha = 0.01) +
#   stat_summary(aes(label=sprintf("%1.1f", ..y..)),
#                geom="text", 
#                fun = function(y) boxplot.stats(y)$stats[c(1,5)],
#                position=position_nudge(x=0.5), 
#                size=3.5) +
#   theme_minimal() +
#   labs(x = "translation frame", y = "jaccard similarity in peptide database")
