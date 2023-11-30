#making supplementary tables for the paper
#Author: Yiqi Huang

library(tidyverse)
supp_table_path <- "/pool/analysis/yhuang/yhuang_iav/reports/supplementary_tables"

#Supplementary Table 1, raw LFQ intensities, simulated intensities and UMAP clusters (+enrichment analysis)----
load("/pool/analysis/yhuang/yhuang_iav/data/phubel_pulse_20171025.RData")
load("/pool/analysis/yhuang/yhuang_iav/results/phubel_pulse_20171025_fit_all.RData")

fasta_info.df <- ms_data$protgroups %>% 
  mutate(gene_names = coalesce(gene_names, str_extract(fasta_headers, "(?<=(GN=))[:graph:]*(?=(\\s|$))")),
         protein_names = coalesce(protein_names, str_extract(fasta_headers, "(?<=(_HUMAN\\s)).*?(?=(\\sOS=))")) )

table1_LFQ.df <- ms_data$protgroup_intensities %>% 
  select(protgroup_id, mschannel, LFQ_intensity = Intensity) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs,gene_names, protein_names, 
                   is_viral, is_contaminant, is_reverse)) %>% 
  left_join(select(ms_data$mschannels, mschannel, mstag, condition, timepoint, replicate)) %>% 
  select(protgroup_id, gene_names, protein_names, majority_protein_acs, is_viral, is_contaminant, is_reverse, 
         condition, replicate, mstag, timepoint, LFQ_intensity) %>% 
  mutate(condition = str_replace(condition, "SC35M", "IAV")) %>% 
  pivot_wider(names_from = timepoint, values_from = LFQ_intensity)

write_tsv(table1_LFQ.df, file = file.path(supp_table_path, "supp_table_1_tab_1_LFQ.txt"))

cluster.df <- read_tsv("/pool/analysis/yhuang/yhuang_iav/data/phubel_pulse_clusters_20190131.txt")
table1_sim_UMAP.df <- fit_stats$label_sim %>% 
  filter(timepoint %in% unique(ms_data$mschannels$timepoint)) %>% 
  mutate(median = `50%`*exp(global_labu_shift)) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs,gene_names, protein_names, 
                   is_viral, is_contaminant, is_reverse)) %>%
  left_join(select(cluster.df, protgroup_id, cluster)) %>% 
  select(protgroup_id, gene_names, protein_names, majority_protein_acs, cluster, is_viral, is_contaminant, is_reverse, 
         condition, mstag, timepoint, median) %>% 
  mutate(condition = str_replace(condition, "SC35M", "IAV")) %>% 
  pivot_wider(names_from = timepoint, values_from = median)

write_tsv(table1_sim_UMAP.df, file = file.path(supp_table_path, "supp_table_1_tab_2_sim_UMAP.txt"))  

#Supplementary Table 2, simulated rates and rate contrasts (+enrichment analysis)----
cum_rates.df <- fit_stats$rate_params %>% 
  filter(var %in% c("s_cum", "d0_cum", "d1_cum", "s_cum_scaled")) %>% 
  select(protgroup_id, condition, rate, param, median = `50%`) %>% 
  pivot_wider(names_from = param, values_from = median) 

table2_sim_rates.df <- fit_stats$rate_sim %>% 
  filter(timepoint %in% unique(ms_data$mschannels$timepoint)) %>% 
  mutate(rate = str_remove_all(var, "_sim")) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs,gene_names, protein_names, 
                   is_viral, is_contaminant, is_reverse)) %>%
  select(protgroup_id, gene_names, protein_names, majority_protein_acs, is_viral, is_contaminant, is_reverse, 
         condition, timepoint, rate, median = `50%`) %>% 
  pivot_wider(names_from = timepoint, values_from = median) %>% 
  left_join(cum_rates.df) %>% 
  mutate(condition = str_replace(condition, "SC35M", "IAV"))

write_tsv(table2_sim_rates.df, file = file.path(supp_table_path, "supp_table_2_tab_1_sim_rates.txt")) 

L0abu.df <- fit_stats$global %>% 
  filter(var == "L0_labu") %>% 
  mutate(median.L0 = pmax(-5, `50%`),
         median.L0.exp = exp(median.L0)) %>% 
  select(protgroup_id, median.L0.exp)

table2_rate_contrasts.df <- select(fit_cond_contrast_stats$rate_params, 
                                   protgroup_id, contrasts, param, rate, difference = `50%`, p_value, sd) %>% 
  filter(!(rate != "s"&str_detect(param, "scale"))) %>%
  left_join(L0abu.df) %>% 
  mutate(difference = ifelse(rate == "s", difference/median.L0.exp, difference),
         sd = ifelse(rate == "s", sd/median.L0.exp, sd)) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs,gene_names, protein_names, 
                   is_viral, is_contaminant, is_reverse)) %>% 
  mutate(contrasts = str_replace_all(contrasts, "SC35M", "IAV"),
         is_sig = p_value <= 0.001,
         is_hit = ifelse(rate == "s", abs(difference) >= 1, abs(difference) >= 0.25)&is_sig&!is_contaminant&!is_reverse,
         rate = ifelse(str_detect(param, "scaled"), paste0(rate, "_scaled"), rate),
         change = case_when(is_hit& difference > 0 ~ "+",
                            is_hit& difference < 0 ~ "-",
                            TRUE ~ ".")) %>% 
  select(protgroup_id, gene_names, protein_names, majority_protein_acs, is_viral, is_contaminant, is_reverse, 
         contrasts, rate, is_hit, change, difference, p_value, sd) %>% 
  pivot_wider(names_from = "rate", values_from = c("is_hit", "change",  "difference", "p_value", "sd"),
               names_sep=".")

write_tsv(table2_rate_contrasts.df, file = file.path(supp_table_path, "supp_table_2_tab_2_rate_contrasts.txt")) 

#Supplementary Table 3, intersection with APMS and genome-wide screen studies----
#apms data
apms_wide.df <- read_tsv(file.path("/pool/analysis/yhuang/yhuang_iav/data/", "metaanalysis", "apms.comb.matched.txt")) %>% 
  replace_na(list(bait = "NA")) %>% 
  mutate(bait = str_replace(bait, "NEP", "NS2"),
         interactor = TRUE) %>% 
  select(-ID) %>% unique() %>% 
  pivot_wider(names_from = screen, values_from = interactor) %>% 
  mutate(across(everything(), ~replace_na(.x, FALSE)),
         total_ident = rowSums(across(heaton:wang)))

turnover4intersection.df <- table2_rate_contrasts.df %>% 
  separate(majority_protein_acs, into = c("protein_ac", "others"), sep = ";") %>% 
  mutate(protein_ac_strip = str_remove_all(protein_ac, "\\-\\d+")) %>% 
  group_by(protein_ac_strip) %>% 
  mutate(sum_detection = n()/3) %>% ungroup() %>% 
  mutate(protein_ac = ifelse(sum_detection >1, protein_ac, protein_ac_strip)) %>% #if the canonical and isoforms are detected separately, keep them separate. Otherwise, always use the canonical ac
  filter(contrasts == "IAV_vs_Mock")

#rematched some proteins inside the protein groups
IAV_contrast.df <- read_tsv("/pool/analysis/yhuang/yhuang_iav/reports/supplementary_tables/supp_table_2_tab_2_rate_contrasts.txt")

IAV_cotrast4intersection.df <- IAV_contrast.df %>% 
  filter(contrasts == "IAV_vs_Mock") %>% 
  separate_rows(majority_protein_acs, sep = ";")

IAV_interactome.df <- read_tsv("/pool/analysis/yhuang/yhuang_iav/reports/supplementary_tables/supp_table_3_tab_1_apms.txt") %>% 
  select(viral_protein:protgroup_id) %>% 
  replace_na(list(viral_protein = "NA"))

IAV_interactome_unmatched.df <- IAV_interactome.df %>% 
  filter(is.na(protgroup_id))

IAVunmatched.df <- IAV_interactome_unmatched.df %>% 
  select(viral_protein, protein_ac, total_ident, gene_name) %>% 
  left_join(IAV_cotrast4intersection.df, by = c("protein_ac" = "majority_protein_acs")) %>% 
  filter(!is.na(protgroup_id)) %>% 
  select(protein_ac, protgroup_id_unmatched = protgroup_id) %>% 
  unique()

IAV_interactome_2.df <- IAV_interactome.df %>% 
  left_join(IAVunmatched.df, by = "protein_ac") %>% 
  mutate(protgroup_id = coalesce(protgroup_id, protgroup_id_unmatched))

IAV_interactome_unmatched_2.df <- IAV_interactome_2.df %>% 
  filter(is.na(protgroup_id)) %>% 
  select(protein_ac, gene_name) %>% 
  unique()

IAV_contrast4intersection_genename.df <- IAV_contrast.df %>% 
  filter(contrasts == "IAV_vs_Mock") %>% 
  separate_rows(gene_names, sep = ";") %>% 
  select(protgroup_id, gene_name = gene_names, majority_protein_acs, contains("change"))

IAVunmatched2.df <- IAV_interactome_unmatched_2.df %>% 
  left_join(IAV_contrast4intersection_genename.df) %>% 
  filter(!is.na(protgroup_id), !str_detect(majority_protein_acs, "REV")) %>% 
  filter(!protgroup_id %in% c(3908, 2582, 2991, 4470, 1508, 2486, 2487, 2489)) %>% #duplicated match, keep the ones with rate changes
  select(protein_ac, protgroup_id_unmatched_2 = protgroup_id)

IAV_interactome_final.df <- IAV_interactome_2.df %>% 
  left_join(IAVunmatched2.df, by = "protein_ac") %>% 
  mutate(protgroup_id = coalesce(protgroup_id, protgroup_id_unmatched_2)) %>% 
  select(-contains("unmatched"))

#reassemble the table 3 apms tab
table3_apms_intersection.df <- IAV_interactome_final.df %>%  
  left_join(select(turnover4intersection.df,  protgroup_id, contains("change"))) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs)) %>% 
  select(viral_protein, protein_ac, gene_name, heaton:majority_protein_acs)

write_tsv(table3_apms_intersection.df, file = file.path(supp_table_path, "supp_table_3_tab_1_apms_new.txt")) 

#genome-wide screens
functional.df <- read_tsv(file.path("/pool/analysis/yhuang/yhuang_iav/data/", "metaanalysis", "all_functional_hits_matched.txt"))

table3_functional_intersection.df <- functional.df %>% 
  mutate(type = ifelse(type == "proviral", 1, -1)) %>% 
  unique() %>% 
  pivot_wider(names_from = source, values_from = type) %>% 
  mutate(score = rowSums(across(Yi:Brass), na.rm = T)) %>% 
  left_join(select(turnover4intersection.df, protein_ac, protgroup_id, contains("change"))) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs, is_contaminant)) %>% 
  select(protein_ac, gene_name, score, Yi:is_contaminant)

functional_unmatched.df <- table3_functional_intersection.df %>% 
  filter(is.na(protgroup_id)) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs) %>% 
              separate_rows(majority_protein_acs, sep = ";") %>% dplyr::rename(protein_ac = majority_protein_acs,
                                                                        protgroup_id_unmatched = protgroup_id))
functional_unmatched_2.df <- filter(functional_unmatched.df, is.na(protgroup_id_unmatched)) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs, gene_names) %>% 
              separate_rows(gene_names, sep = ";") %>% rename(gene_name = gene_names,
                                                                        protgroup_id_unmatched_2 = protgroup_id))
#There was no new matching by gene name
table3_functional_intersection_final.df <- table3_functional_intersection.df %>% 
  left_join(functional_unmatched.df %>% select(protein_ac, protgroup_id_unmatched)) %>% 
  mutate(protgroup_id = coalesce(protgroup_id, protgroup_id_unmatched)) %>% 
  select(protein_ac:protgroup_id) %>% 
  left_join(select(turnover4intersection.df, protgroup_id, contains("change"))) %>% 
  left_join(select(fasta_info.df, protgroup_id, majority_protein_acs, is_contaminant)) %>% 
  select(protein_ac, gene_name, score, Yi:is_contaminant)

write_tsv(table3_functional_intersection_final.df, file = file.path(supp_table_path, "supp_table_3_tab_2_functional_new.txt"))

#Supplementary Table 4, raw data from functional screens----
#reused some tables made during figure generation
table4_renilla.df <- kd_screen_72.df %>% 
  filter((cell.line == "hela"&biorep == 2)|
           (cell.line == "a549"&gene_name %in% c("GPKOW", "JAK1", "STAT3", "CLIC4", "Scrambled"))) %>% 
  filter(!((gene_name == "LGALS3BP")&plate != 3)) %>% 
  left_join(select(kd_screen_72_ttest.df,
                   cell.line, plate, gene_name = group2, biorep,log2FC = estimate, p, p.adj, p.adj.signif)) %>% 
  left_join(turnover.df %>% filter((turnover_category != "."&contrasts != "SC35MdelNS1_vs_Mock")|gene_name == "JAK1") %>%
              select(protgroup_id, gene_name) %>% unique()) %>% 
  mutate(log2FC = -log2FC,
         type = case_when(is_control ~ "control",
                          p.adj.signif != "ns" & log2FC > 0 ~ "antiviral",
                          p.adj.signif != "ns" & log2FC < 0 ~ "proviral",
                          TRUE ~ "not significant"),
         type = factor(type, levels = c("control", "antiviral", "proviral", "not significant"))) %>% 
  select(cell_line = cell.line, virus_strain = virus, metric, 
         gene_name, protgroup_id, replicate = well_number, plate, intensity, intensity_log2,
         log2FC, p_value = p, p_value_adj = p.adj, significance = p.adj.signif, type) %>% 
  mutate(cell_line = str_replace(cell_line, "hela", "HeLa"),
         cell_line = str_replace(cell_line, "a549", "A549"),
         virus_strain = toupper(virus_strain),
         metric = "Renilla") %>% 
  arrange(cell_line, gene_name)

write_tsv(table4_renilla.df, file = file.path(supp_table_path, "supp_table_4_tab_2_renilla.txt"))

table4_mini_replicon.df <- kd_screen_80_relevant.df %>% 
  left_join(kd_screen_80_relevant.df %>% filter(gene_name == "Scrambled") %>% select(cell.line, biorep, plate, ratio_control = ratio_mean) %>% unique()) %>% 
  mutate(ratio_norm = ratio_mean/ratio_control) %>% 
  mutate(ratio_norm_log2 = log2(ratio_norm)) %>% 
  filter(!str_detect(gene_name, "no siRNA")) %>% 
  filter((cell.line == "hela")|
           (cell.line == "A549"&gene_name %in% c("GPKOW", "JAK1", "STAT3", "CLIC4", "Scrambled"))) %>% 
  mutate(effect = case_when(ratio_norm_log2 > 0.5 ~ "+",
                            ratio_norm_log2 < -0.5 ~ "-",
                            TRUE ~ "."),
         intensity_ratio_log2 = log2(ratio_mean)) %>% 
  left_join(select(kd_screen_80_norm_2_wide.df,
                   gene_name, cell.line, type)) %>% 
  left_join(turnover.df %>% filter((turnover_category != "."&contrasts != "SC35MdelNS1_vs_Mock")|gene_name == "JAK1") %>%
              select(protgroup_id, gene_name) %>% unique()) %>% 
  select(cell_line = cell.line, virus_strain = virus, metric, 
         gene_name, protgroup_id, replicate = biorep, plate,  intensity_ratio = ratio_mean, intensity_ratio_log2,
         log2FC = ratio_norm_log2, effect, type) %>% 
  mutate(cell_line = str_replace(cell_line, "hela", "HeLa"),
         #cell_line = str_replace(cell_line, "a549", "A549"),
         virus_strain = toupper(virus_strain),
         metric = "Firefly/Renilla") %>% 
  arrange(cell_line, gene_name)

write_tsv(table4_mini_replicon.df, file = file.path(supp_table_path, "supp_table_4_tab_3_mini_replicon.txt"))

table4_viability.df <- kd_viability_plot.df <- kd_screen_viability4paper.df %>% 
  filter((cell.line == "hela" )|(cell.line == "a549"&
         gene_name %in% c("GPKOW", "JAK1", "STAT3", "CLIC4", "Scrambled"))) %>% 
  left_join(kd_screen_viability4paper_summary.df) %>% 
  left_join(turnover.df %>% filter((turnover_category != "."&contrasts != "SC35MdelNS1_vs_Mock")|gene_name == "JAK1") %>%
              select(protgroup_id, gene_name) %>% unique()) %>% 
  left_join(kd_viability_ttest.df %>% select(cell.line,  plate, gene_name = group2, p_value = p, p_value_adj = p.adj)) %>% 
  select(cell_line = cell.line,  metric, 
         gene_name, protgroup_id, replicate, plate,  
         intensity, #control_intensity_mean, 
         mean_viability = mean, p_value, p_value_adj, type = group) %>% 
  mutate(cell_line = str_replace(cell_line, "hela", "HeLa"),
         cell_line = str_replace(cell_line, "a549", "A549"),
         metric = "Resazurin") %>% 
  arrange(cell_line, gene_name)

write_tsv(table4_viability.df, file = file.path(supp_table_path, "supp_table_4_tab_1_viability.txt"))

#Supplementary Table 5, proteomic analysis of GPKOW KD/KO cells----
response_to_virus.df <- read_tsv(file.path("/pool/analysis/yhuang/yhuang_iav/data/", "pubdata", "response_to_virus.tsv"))
defense_response.df <- response_to_virus.df %>% 
  filter(str_detect(`GO NAME`, "defense")) %>% 
  select(protein_ac = `GENE PRODUCT ID`, gene_name_db = SYMBOL, gobp = `GO NAME`) %>% 
  unique()
load(file.path("/pool/analysis/yhuang/yhuang_iav/results/", "yhuang_iav_msglm_fit_phubel_exp76_20231011b.RData"))

table5_kd_infection.df <- object_contrasts.df %>% 
  filter(ci_target == "average") %>% 
  mutate(protein_ac = str_remove_all(protac_label, fixed("..."))) %>% 
  left_join(defense_response.df) %>% 
  mutate(is_defense_response = !is.na(gobp),
         change = ifelse(is_hit, ifelse(median > 0 , "+", "-"), "."),
         contrasts = "siGPKOW_vs_control@IAV") %>% 
  left_join(fit_stats$objects %>% select(object_id, gene_names, majority_protein_acs, is_viral, is_contaminant, is_reverse) %>% unique()) %>% 
  select(protgroup_id = object_id, gene_names, majority_protein_acs, is_defense_response, is_viral, is_contaminant, is_reverse,
         contrasts, is_hit, change, fold_change_log2 = median, p_value, sd_log2 = sd)

write_tsv(table5_kd_infection.df, file = file.path(supp_table_path, "supp_table_5_tab_1_kd_infection.txt"))

load(file.path("/pool/analysis/yhuang/yhuang_iav/results/", "yhuang_iav_msglm_fit_phubel_exp123_20231011b.RData"))

table5_ntc_treatment.df <- object_contrasts.df %>% 
  filter(ci_target == "average", str_detect(contrast, "_vs_Mock@NTC")) %>% 
 mutate(protein_ac = str_remove_all(protac_label, fixed("...")),
         contrast = str_remove_all(contrast, "@NTC")) %>% 
  left_join(defense_response.df) %>% 
  mutate(is_defense_response = !is.na(gobp),
         change = ifelse(is_hit, ifelse(median > 0 , "+", "-"), ".")) %>% 
  left_join(fit_stats$objects %>% select(object_id, gene_names, majority_protein_acs, is_viral, is_contaminant, is_reverse) %>% unique()) %>% 
  select(protgroup_id = object_id, gene_names, majority_protein_acs, is_defense_response, is_viral, is_contaminant, is_reverse,
         contrast, is_hit, change, fold_change_log2 = median, p_value, sd_log2 = sd) %>% 
  pivot_wider(names_from = "contrast", values_from = c("is_hit", "change",  "fold_change_log2", "p_value", "sd_log2"),
              names_sep=".")

write_tsv(table5_ntc_treatment.df, file = file.path(supp_table_path, "supp_table_5_tab_2_ntc_treatment.txt"))

table5_koXtreatment.df <-  object_contrasts.df %>% 
  filter(ci_target == "average", str_detect(contrast, "gpkow_vs_ntc")) %>% 
  mutate(protein_ac = str_remove_all(protac_label, fixed("...")),
         contrast = str_remove_all(contrast, "gpkow_vs_ntc@")) %>% 
  left_join(defense_response.df) %>% 
  mutate(is_defense_response = !is.na(gobp),
         change = ifelse(is_hit, ifelse(median > 0 , "+", "-"), ".")) %>% 
  left_join(fit_stats$objects %>% select(object_id, gene_names, majority_protein_acs, is_viral, is_contaminant, is_reverse) %>% unique()) %>% 
  select(protgroup_id = object_id, gene_names, majority_protein_acs, is_defense_response, is_viral, is_contaminant, is_reverse,
         contrast, is_hit, change, fold_change_log2 = median, p_value, sd_log2 = sd) %>% 
  pivot_wider(names_from = "contrast", values_from = c("is_hit", "change",  "fold_change_log2", "p_value", "sd_log2"),
              names_sep=".")

write_tsv(table5_koXtreatment.df, file = file.path(supp_table_path, "supp_table_5_tab_3_koXtreatment.txt"))