library(dplyr)
library(readr)
library(tidyr)
library(taxize)
setwd("github/2021-orpheum-sim/")

# find refseq not in GTDB -------------------------------------------------
gtdb <- read_csv("https://osf.io/p6z3w/download")
gtdb <- gtdb %>%
  mutate(accession_minus_version = gsub("\\.*", "", ident)) %>%
  mutate(accession_minus_prefix  = gsub("GC[FA]_", "", accession_minus_version))

refseq <- read_tsv("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt",
                   skip = 1) 
refseq <- refseq %>%
  rename(assembly_accession = `# assembly_accession`) %>%
  mutate(accession_minus_version = gsub("\\.*", "", assembly_accession)) %>%
  mutate(accession_minus_prefix  = gsub("GC[FA]_", "", accession_minus_version))

new_refseq <- refseq %>%
  filter(! accession_minus_prefix %in% gtdb$accession_minus_prefix)  %>%
  filter(seq_rel_date > "2021-04-27") %>% # filter to genomes released after GTDB rs202
  filter(!is.na(infraspecific_name))      # filter genomes with no infraspecific name; most are phages or euks

tmp <- new_refseq %>%
  select(organism_name) %>%
  group_by(organism_name) %>%
  tally()

# try and get lineages to chose subset of new genomes -----------------------

uid_to_lin <- classification(sci_id = new_refseq$taxid, db = "ncbi", 
                             return_id = T, batch_size = 50, max_tries = 2)
uid_to_lin_df <- data.frame(superkingdom = character(),
                            phylum = character(), 
                            class  = character(),
                            order  = character(),
                            family = character(),
                            genus  = character(), 
                            species= character(),
                            id     = character())
for(i in 1:length(uid_to_lin)){
  uid_to_lin_n <- uid_to_lin[[i]] %>%
    select(name, rank) %>%
    filter(rank %in% c("superkingdom", "phylum","class","order","family", "genus", "species")) %>%
    pivot_wider(names_from = rank, values_from = name) %>%
    mutate(id = names(uid_to_lin[i]))
  uid_to_lin_df <- bind_rows(uid_to_lin_df, uid_to_lin_n)
}

uid_to_lin_df <- uid_to_lin_df %>%
  distinct() 

# write_tsv(uid_to_lin_df, "new_refseq_lineages.tsv")

# select new refseq representatives ---------------------------------------

new_refseq <- new_refseq %>%
  mutate(taxid = as.character(taxid)) %>%
  left_join(uid_to_lin_df, by = c("taxid" = "id")) %>%
  filter(superkingdom %in% c("Archaea", "Bacteria")) 

new_refseq %>% 
  group_by(superkingdom, phylum) %>%
  tally()

new_refseq_slice <- new_refseq %>%
  group_by(phylum) %>%
  slice(1) %>%
  mutate(species_no_space = gsub(" ", "-", species))

write_tsv(new_refseq_slice, "inputs/refseq_not_in_gtdb_metadata.tsv")

# select representatives in GTDB ------------------------------------------

# select the same set of species used in 2021-panmers (only non-IBD species)

in_gtdb <- gtdb %>%
  group_by(superkingdom, phylum, class, order, family, genus, species) %>%
  tally() %>%
  filter(n >= 20) %>%
  filter(n < 1000) %>% # temporarily remove very large pangenomes for testing.
  arrange(desc(n))

in_gtdb_slice <- in_gtdb %>%
  group_by(phylum) %>%
  slice(1) %>%
  mutate(species_no_space = gsub(" ", "-", species))

# determine representative _genome_ in GTDB from set of speices

gtdb_metadata_bacteria <- read_tsv("https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_metadata_r202.tar.gz") %>%
  rename(ident = bac120_metadata_r202.tsv)
gtdb_metadata_archaea <- read_tsv("https://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_metadata_r202.tar.gz") %>%
  rename(ident=ar122_metadata_r202.tsv)

gtdb_rep_metadata <- bind_rows(gtdb_metadata_bacteria, gtdb_metadata_archaea) %>%
  filter(gtdb_representative == TRUE)

gtdb_rep_metadata_slice <- gtdb_rep_metadata %>%
  separate(gtdb_taxonomy, into = c('superkingdom', 'phylum', 'class', 'order', 
                                   'family', 'genus', 'species'), sep = ";") %>%
  filter(species %in% in_gtdb_slice$species) %>%
  mutate(ident = gsub("GB_", "", ident)) %>%
  mutate(ident = gsub("RS_", "", ident)) %>%
  mutate(species_no_space = gsub(" ", "-", species))
  
write_tsv(gtdb_rep_metadata_slice, "inputs/gtdb_metadata.tsv")


# combine not in GTDB and GTDB into one metadata file ---------------------

gtdb_rep_metadata_slice2 <- gtdb_rep_metadata_slice %>%
  mutate(accession_minus_version = gsub("\\.*", "", ident)) %>%
  mutate(accession_minus_prefix  = gsub("GC[AF]_", "", accession_minus_version))


genbank <- read_tsv("https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt",
                    skip = 1)
genbank <- genbank %>%
  mutate(accession_minus_version = gsub("\\.*", "", `# assembly_accession`)) %>%
  mutate(accession_minus_prefix  = gsub("GC[AF]_", "", accession_minus_version))



# add the extra accession that's just not in the spreadsheet but is in genbank...
extra <- tibble(`# assembly_accession`  = "GCA_000635495.1",
                biosample               = "SAMN02671343",
                bioproject              = "PRJNA239872",
                taxid                   = 1471472,
                organism_name           = "Prochlorococcus sp. scB243_495K23 (cyanobacteria)",
                ftp_path = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/635/495/GCA_000635495.1_De_novo_assembly_of_single-cell_genomes/',
                accession_minus_version = "GCA_000635495",
                accession_minus_prefix  = "000635495")

gtdb_genbank <- genbank %>%
  filter(accession_minus_prefix %in% gtdb_rep_metadata_slice2$accession_minus_prefix) %>%
  bind_rows(extra) %>%
  mutate(set = "gtdb_representatives")

refseq_genbank <- genbank %>%
  filter(accession_minus_prefix %in% new_refseq_slice$accession_minus_prefix) %>%
  mutate(set = "refseq_not_in_genbank")
  
all_genomes <- bind_rows(gtdb_genbank, refseq_genbank) 
colnames(all_genomes) <- c('assembly_accession', 'bioproject', 'biosample',
                           'wgs_master', 'refseq_category', 'taxid', 'species_taxid', 
                           'organism_name', 'infraspecific_name', 'isolate', 
                           'version_status', 'assembly_level', 'release_type', 
                           'genome_rep', 'seq_rel_date', 'asm_name', 'submitter', 
                           'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path', 
                           'excluded_from_refseq', 'relation_to_type_material', 
                           'asm_not_live_date', 'accession_minus_version', 
                           'accession_minus_prefix', 'set')

write_tsv(all_genomes, "inputs/all_genomes_genbank_info_metadata.tsv")

# name outputs by accession, deal with taxonomy/species names later;
# can be rescued by left_join() with the other metadata tables writeen above
