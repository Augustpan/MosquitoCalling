library(tidyverse)

names_blast = c("lib_id", "contig", "hit", "pident", "alnlen", "evalue", "bitscore",
                "phylum", "class", "order", "family", "genus", "species")
tophits = read_tsv("intermediate/all_contigs_tophits.txt", col_names=names_blast)

coi = read_lines("intermediate/coi_contigs.list")
tophits = filter(tophits, contig %in% coi)

names_snp = c("lib_id", "snp_count", "mean_frequency", "mean_depth")
snp = read_tsv("intermediate/stat_snp.txt", col_names=names_snp)

lib_ids = unique(filter(tophits, family=="Culicidae")$lib_id)
tophits_mos = filter(tophits, family=="Culicidae")
identified = NULL
contaminated = NULL
todo = NULL
for (lid in lib_ids) {
    x = filter(tophits_mos, lib_id==lid)
    sp_list = unique(x$species) %>% sort()
    if (any(x$pident < 98)) {
        if (length(sp_list) == 1) {
            todo = rbind(todo, c(lid, paste0(sp_list, collapse = "/"), "no_but_species_uncertain"))
        } else {
            todo = rbind(todo, c(lid, paste0(sp_list, collapse = "/"), "unsure"))
        }
    } else {
        if (length(sp_list) > 1) {
            contaminated = rbind(contaminated, c(lid, paste0(sp_list, collapse = "/"), "yes"))
        } else {
            identified = rbind(identified, c(lid, sp_list, "no"))
        }
    }
}

all = rbind(todo, contaminated, identified)

colnames(all) = c("lib_id", "species", "is_contaminated")
all = as_tibble(all) %>% left_join(snp)
all_merged = filter(tophits, family!="Culicidae") %>% 
    distinct(lib_id) %>% 
    mutate(contains_non_Culicidae_COI="yes") %>%
    right_join(all) %>%
    replace_na(list(other_coi="no"))

filter(all_merged, is_contaminated=="unsure")$lib_id %>% 
    write_lines("intermediate/uncertain_libs.txt")
write_csv(all_merged, "intermediate/all.csv", na="")
write_lines(tophits_mos$contig, "intermediate/culicidae_contigs.txt")
