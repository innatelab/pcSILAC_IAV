# Prepare data for the full proteome analysis for Philipp's exp123 - GPKOW KO A549 cells + IVT4 (IFNbd doesn't have mock control)
# Experiments done by Philipp Hubel in 2017
# Author: Yiqi Huang
###############################################################################

project_id <- 'yhuang_iav'
message('Project ID=', project_id)
data_version <- "20231010"
fit_version <- "20231011b" #b uses the more stringent normalisation 
mstype <- "phubel_exp123"
message('Dataset version is ', data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, "ggplot_ext.R"))

library(msimportr)
library(msglm)
library(tidyverse)
library(jsonlite)
library(pheatmap)

msfolder <- "phubel_exp76+123" #this was quanted together with exp76
msdata_path <- file.path(data_path, str_c(msfolder,"_", data_version))

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  pepmodstate_mscalib_filename = "mscalib_QEP5_intensity_pepmodstate_cov2_20211108.json",
                  msfolder = msfolder, quant_type = "LFQ", quant_col_prefix = "LFQ_Intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$pepmodstate_mscalib_filename, '...')
pepmodstate_mscalib <- read_mscalib_json(file.path(data_path, data_info$pepmodstate_mscalib_filename)) 

msruns.df <- read_tsv(file.path(msdata_path, "combined/experimentalDesign.txt")) %>% 
  dplyr::rename(rawfile=Name, fraction=Fraction, msexperiment_mq=Experiment, is_ptm = PTM) %>%
  tidyr::separate("msexperiment_mq", into = c("cell_type", "ko", "treatment", "replicate"),
                 sep = "_", remove = FALSE) %>% 
  dplyr::mutate(ko = factor(ko, c("NTC", "GPKOW", "Scrambled", "siGPKOW")),
                treatment = factor(treatment, c("mock", "IVT4", "IFN", "IAV")),
                condition = str_c(cell_type, "_", ko,"_", treatment),
                replicate = as.integer(replicate),
                is_skipped = cell_type == "HeLa") %>% 
  dplyr::arrange(ko, treatment, replicate, cell_type) %>% 
  dplyr::mutate(condition = factor(condition, levels=unique(condition)),
                msexperiment = str_c(condition, "_", replicate),
                msexperiment = factor(msexperiment, levels=unique(msexperiment))
  ) %>% 
  dplyr::arrange(ko, treatment, replicate)

fasta.dfs <- list(
  human = read_innate_uniprot_fasta(file.path(msdata_path, "fasta/2023_04_uniprotkb_proteome_UP000005640_AND_revi_2023_10_05.fasta")),
  IAV = read_innate_uniprot_fasta(file.path(msdata_path, "fasta/2023_04_uniprotkb_proteome_UP000008576_2023_10_05.fasta"))
)

msdata.wide <- read.MaxQuant.ProteinGroups(file.path(msdata_path, 'combined/txt'), import_data = c(data_info$quant_type, "ident_type"))
msdata_colgroups <- attr(msdata.wide, "column_groups")

mqevidence <- read.MaxQuant.Evidence(file.path(msdata_path, 'combined', 'txt'),
                                     mschannel_annotate.f = function(mschans_df) {
                                       res <- dplyr::inner_join(dplyr::select(mschans_df, mstag, rawfile),
                                                                msruns.df,
                                                                by="rawfile")
                                       attr(res, "column_scopes") <- c(ko = "msexperiment",
                                                                       treatment = "msexperiment",
                                                                       condition = "msexperiment",
                                                                       replicate = "msexperiment")
                                       return(res)
                                     })
mqevidence$peptides <- read.MaxQuant.Peptides(file.path(msdata_path, 'combined', 'txt'), file_name='peptides.txt',
                                              import_data='ident_type')
mqevidence$peaks <- NULL # exclude big data frame

strlist_label <- function(strs) {
  str_c(strs[[1]], if_else(n_distinct(strs) > 1, '...', ''))
}
strlist_label2 <- function(strs, delim=fixed(';')) {
  sapply(str_split(strs, delim), strlist_label)
}

# all ms data
msdata_full <- list(msexperiments = dplyr::select(mqevidence$msexperiments, -starts_with("msfraction")),
                    msruns = mqevidence$msruns)

msdata_full <- append_protgroups_info(msdata_full, dplyr::mutate(msdata.wide, organism=NULL),
                                      proteins_info = dplyr::bind_rows(fasta.dfs) %>%
                                        dplyr::mutate(protein_ac_noiso = str_remove(protein_ac, "-\\d+(?:#.+)*$"),
                                                      protein_isoform_ix = replace_na(as.integer(str_match(protein_ac, "-(\\d+)$")[, 2]), 1L)),
                                      import_columns = c("organism"))


msdata_full$peptides <- mqevidence$peptides %>%
  dplyr::left_join(select(dplyr::bind_rows(fasta.dfs), lead_razor_protein_ac = protein_ac, organism)) %>%
  dplyr::mutate(peptide_rank = 1L)
msdata_full$pepmods <- mqevidence$pepmods %>%
  dplyr::mutate(pepmod_rank = 1L)
msdata_full$pepmodstates <- mqevidence$pepmodstates
# pepmods and peptides could be ranked by their relevance to the experiment, e.g. APMS, compartment... see examples in adhoc 

# redefine protein groups (protregroups)
peptides.df <- dplyr::select(msdata_full$peptides, peptide_id, protgroup_ids, protein_acs, lead_razor_protein_ac, peptide_seq, is_reverse, peptide_rank)
proteins.df <- msdata_full$proteins
save(file = file.path(msdata_path, str_c(project_id, "_", msfolder, "_peptides.RData")),
     peptides.df, proteins.df)

# .. run protregroup.jl
msdata_full$protregroups <- read_tsv(file.path(msdata_path, 
                                               str_c(project_id, "_", msfolder,"_", data_version,  "_protregroups.txt")),
                                     col_types = list(protregroup_id = "i"))

msdata_full$protein2protregroup <- dplyr::select(msdata_full$protregroups, protregroup_id, protein_ac=majority_protein_acs) %>%
  separate_rows(protein_ac, sep=fixed(";"), convert=TRUE) %>%
  dplyr::mutate(is_majority = TRUE) %>%
  dplyr::group_by(protregroup_id) %>%
  dplyr::mutate(protein_ac_rank = row_number()) %>%
  dplyr::ungroup()

msdata_full$protregroup2peptide <- bind_rows(
  select(msdata_full$protregroups, protregroup_id, peptide_id=spec_peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = TRUE),
  select(msdata_full$protregroups, protregroup_id, peptide_id=peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = FALSE)) %>%
  dplyr::group_by(protregroup_id, peptide_id) %>%
  dplyr::summarise(is_specific = any(is_specific)) %>%
  dplyr::ungroup()
msdata_full$protregroup2pepmod <- dplyr::inner_join(msdata_full$protregroup2peptide,
                                                    dplyr::select(msdata_full$pepmods, pepmod_id, peptide_id)) %>%
  dplyr::select(-peptide_id)

msdata_full$protregroup2pepmodstate <- dplyr::semi_join(msdata_full$protregroup2pepmod,
                                                        dplyr::select(msdata_full$pepmods, pepmod_id)) %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmodstates, pepmod_id, pepmodstate_id))

msdata_full$protregroups <- msdata_full$protregroups %>%
  dplyr::mutate(gene_label = strlist_label2(gene_names),
                protein_label = strlist_label2(protein_names),
                protein_description = strlist_label2(protein_descriptions),
                is_viral = replace_na(str_detect(organism, "virus"), FALSE),
                protac_label = strlist_label2(majority_protein_acs),
                protregroup_label = case_when(is_viral ~ protein_label,
                                              !is.na(gene_label) ~ gene_label,
                                              !is.na(protac_label) ~ protac_label,
                                              TRUE ~ str_c('#', protregroup_id))) %>%
  dplyr::left_join(dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmods) %>%
                     dplyr::group_by(protregroup_id) %>%
                     dplyr::summarise(npeptides = n_distinct(peptide_id),
                                      npepmods = n_distinct(pepmod_id),
                                      npeptides_unique = n_distinct(peptide_id[is_specific]),
                                      npepmods_unique = n_distinct(pepmod_id[is_specific])) %>%
                     dplyr::ungroup() %>%
                     dplyr::mutate(npeptides_unique_razor = npeptides_unique,
                                   npeptides_razor = 0L,
                                   npepmods_unique_razor = npepmods_unique,
                                   npepmods_razor = 0L))

msdata_full$pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::mutate(is_idented = str_detect(ident_type, "MSMS"))

# prepare protgroup intensities (wider format: all mstags in one row)
intensity_prespec_df <- tibble(.name = msdata_colgroups$LFQ) %>%
  extract(.name, c("mstag", "msexperiment"), remove=FALSE,
          str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)")) %>%
  mutate(.value = str_c("intensity.", mstag)) %>% 
  dplyr::inner_join(select(msruns.df, msexperiment, rawfile))

# prepare protgroup intensities (longer format: each mstag on its own row)
protgroup_intensities_all.df <- tidyr::pivot_longer_spec(
  dplyr::select(msdata.wide, protgroup_id, !!msdata_colgroups$LFQ),
  mutate(intensity_prespec_df, .value = "intensity")) %>%
  select(-mstag)

msdata_full$protgroup_intensities <- protgroup_intensities_all.df %>%
  dplyr::semi_join(select(msdata_full$msruns, msexperiment, rawfile)) %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::select(-rawfile)

msdata_full$protregroup_idents <- dplyr::inner_join(msdata_full$protregroup2pepmod,
                                                    msdata_full$pepmodstate_intensities) %>%
  dplyr::group_by(msexperiment, protregroup_id) %>%
  dplyr::summarise(npepmods_quanted = sum(!is.na(intensity)),
                   nspecpepmods_quanted = sum(!is.na(intensity) & is_specific),
                   npepmods_idented = sum(is_idented),
                   nspecpepmods_idented = sum(is_idented & is_specific)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ident_type = factor(if_else(nspecpepmods_idented > 0, "By MS/MS", "By matching")))

# set up experimental design matrices
# condition = koXtreatment
conditions.df <- dplyr::select(msdata_full$msruns, condition, ko, treatment) %>%
  droplevels() %>% 
  dplyr::distinct() 

conditionXeffect.mtx <- model.matrix(~1 + ko*treatment, data = conditions.df)
conditionXeffect.mtx <- conditionXeffect.mtx[, (apply(conditionXeffect.mtx, 2, function(x) min(abs(x))) == 0)]
dimnames(conditionXeffect.mtx) <- list(condition = as.character(conditions.df$condition),
                                       effect = colnames(conditionXeffect.mtx))
all_conditions <- as.character(conditions.df$condition)

plots_path <- file.path(analysis_path, "plots", str_c(mstype, "_",data_version,"_", fit_version))
if (!dir.exists(plots_path)) dir.create(plots_path)

pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(plots_path,
                              paste0(project_id, "_", msfolder, "_exp_design_",fit_version, ".pdf")),
         width = 6, height = 6)

effects.df <- tibble(effect = colnames(conditionXeffect.mtx)) %>% 
  mutate(treatment = effect_factor(effect, "treatment", levels(conditions.df$treatment), NA),
         ko = effect_factor(effect, "ko", levels(conditions.df$ko), NA),
         effect_type = case_when(
           !is.na(ko) & !is.na(treatment) ~ "ko_treatment",
           !is.na(ko) ~ "ko",
           !is.na(treatment) ~ "treatment",
           TRUE ~ NA_character_),
         prior_mean = 0, #We don't expect any general trends due to any treatment, so the mean is 0 here.
         prior_tau = case_when(effect_type == "treatment" ~ 1.0,
                               effect_type == "ko" ~ 1.0,
                               effect_type == "ko_treatment" ~ 0.25,
                               TRUE ~ 1.0), #This part belongs to the horseshoe prior. Tau represents expected none-zero values. Controls how much the model would shrink the "low values". The higher it is, the more none-zero effect we expect and the less it will shrink!!! 
         prior_df1 = case_when(effect_type == "treatment" ~ 2.0,
                               effect_type == "ko" ~ 2.0,
                               effect_type == "ko_treatment" ~ 4.0,
                               TRUE ~ 1.0), #This part controls the horseshoe+ prior. Controls how conservative the model is. We can leave everything at 2.0 (default value) or just delete it, and adjust it later if we see problems in the outcome.
         prior_df2 = prior_df1, 
         is_positive = FALSE)

effects.df$effect_label <- sapply(1:nrow(effects.df), function(i) {
  ko <- as.character(effects.df$ko[[i]])
  trt <- as.character(effects.df$treatment[[i]])
  chunks <- if (!is.na(ko)) {
    ko
  } else {
    character(0)
  }
  if (!is.na(trt)) chunks <- append(chunks, trt)
  paste0(chunks, collapse="+")
})

#The metaconditions in this case correspond to the conditions
all_metaconditions <- c(all_conditions)
conditionXmetacondition.mtx <- msglm::constant_matrix(FALSE, list(condition = all_conditions,
                                                                  metacondition = all_metaconditions))

for (cname in levels(conditions.df$condition)) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
}

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>% 
  filter(n != 0) %>% 
  select(-n)

pheatmap(ifelse(conditionXmetacondition.mtx, 1L, 0L), cluster_rows=FALSE, cluster_cols=FALSE,
          filename = file.path(plots_path,
                               paste0(project_id,"_", msfolder, "_metaconditions_",fit_version, ".pdf")),
          width = 8, height = 6)

gpkow_vs_ntc <- paste0("gpkow_vs_ntc@", levels(conditions.df$treatment))
IVT4_vs_mock <- paste0("IVT4_vs_Mock@", levels(conditions.df$ko))
IFN_vs_mock <- paste0("IFN_vs_Mock@", levels(conditions.df$ko))
IVT4_vs_IFN <- paste0("IVT4_vs_IFN@", levels(conditions.df$ko))
all_contrasts <- c(gpkow_vs_ntc, IVT4_vs_mock, IFN_vs_mock, IVT4_vs_IFN)
metaconditionXcontrast.mtx <- msglm::constant_matrix(0, list(metacondition = all_metaconditions,
                                                             contrast = all_contrasts))

for (ko in levels(conditions.df$ko)) {
  metaconditionXcontrast.mtx[paste0('A549_', ko, '_IVT4'), paste0('IVT4_vs_Mock@', ko)] <- 1.0
  metaconditionXcontrast.mtx[paste0('A549_', ko, '_mock'), paste0('IVT4_vs_Mock@', ko)] <- -1.0
}
for (ko in levels(conditions.df$ko)) {
  metaconditionXcontrast.mtx[paste0('A549_', ko, '_IFN'), paste0('IFN_vs_Mock@', ko)] <- 1.0
  metaconditionXcontrast.mtx[paste0('A549_', ko, '_mock'), paste0('IFN_vs_Mock@', ko)] <- -1.0
}
for (ko in levels(conditions.df$ko)) {
  metaconditionXcontrast.mtx[paste0('A549_', ko, '_IVT4'), paste0('IVT4_vs_IFN@', ko)] <- 1.0
  metaconditionXcontrast.mtx[paste0('A549_', ko, '_IFN'), paste0('IVT4_vs_IFN@', ko)] <- -1.0
}
for (trt in levels(conditions.df$treatment)) {
  metaconditionXcontrast.mtx[paste0('A549_GPKOW_', trt), paste0('gpkow_vs_ntc@', trt)] <- 1.0
  metaconditionXcontrast.mtx[paste0('A549_NTC_', trt), paste0('gpkow_vs_ntc@', trt)] <- -1.0
}

pheatmap(metaconditionXcontrast.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plots_path,
                              paste0(project_id,"_", msfolder, "_contrasts_", data_info$fit_ver, ".pdf")),
         width = 6, height = 5)


metaconditionXcontrast.df <- as.data.frame.table(metaconditionXcontrast.mtx, responseName="weight",
                                                 stringsAsFactors = FALSE) %>% as_tibble() %>%
  filter(weight != 0) %>% 
  mutate(contrast_type = "comparison",
         condition_role = if_else(contrast_type == 'filtering',
                                  if_else(weight > 0, 'signal', 'background'),
                                  'signal'))

contrasts.df <- contrasts.df <- dplyr::select(metaconditionXcontrast.df, contrast, contrast_type) %>%
  dplyr::distinct()

#no batch effects, everything done at once

msglm_def <- msglm_model(conditionXeffect.mtx, conditions.df, effects.df,
                         verbose=TRUE) %>%
  msglm::set_contrasts(metaconditionXcontrast.mtx, conditionXmetacondition.mtx,
                       contrasts.df)

msdata <- import_msglm_data(msdata_full, msglm_def, 
                            object="protregroup", quantobject="pepmodstate", verbose = TRUE)

##########
## normalization
msdata4norm.df <- msdata_full$pepmodstate_intensities %>%
  dplyr::semi_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id)) %>%
  dplyr::semi_join(dplyr::semi_join(msdata_full$pepmods,
                                    dplyr::filter(msdata_full$protregroups, !is_reverse & !is_contaminant & !is_viral) %>%
                                      dplyr::inner_join(msdata_full$protregroup2pepmod) %>% dplyr::select(pepmod_id)) %>%
                     dplyr::filter(!is_reverse & !is_contaminant) %>%
                     dplyr::select(pepmod_id))


require(cmdstanr)
options(mc.cores=8L)
Sys.setenv(MKL_NUM_THREADS=1L)

# normalise msruns within the group of msruns with the same condition and msfraction
msruns_hnorm <- multilevel_normalize_experiments(pepmodstate_mscalib,
                                                 msdata$msruns,
                                                 msdata4norm.df,
                                                 quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "msrun",
                                                 #mcmc.iter = 2000L,
                                                 verbose=TRUE,
                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=2000L),
                                                                    condition = list(cond_col="condition", max_objs=1000L),
                                                                    treatment = list(cond_col = "treatment", max_objs=1000L))) #Set the levels so the group with the smallest difference comes first and the largest difference last.
msdata$msrun_shifts <- msruns_hnorm$mschannel_shifts

#boxplot to check normalisation
boxplotdata <-left_join(msdata$pepmodstate_intensities, select(msdata$msruns, msrun)) %>%
  left_join(msdata$msrun_shifts) %>%
  mutate(norm_intensity = intensity*exp(-total_msrun_shift),
         log2intensity = log2(norm_intensity)) %>%
  select(log2intensity, msrun)

(p <- ggplot(boxplotdata, aes(x=msrun, y=log2intensity)) +
    geom_boxplot() + theme(axis.text.x = element_text(angle = 90)))

### more stringent normalisation (used in the "b" versions)
row_median <- msdata4norm.df %>%  
  group_by(pepmodstate_id) %>% summarise(protgroup_median=median(log(intensity))) %>% ungroup()

intensity_row_median <- msdata4norm.df %>% 
  left_join(row_median) %>% mutate(norm_factor_row=log(intensity) - protgroup_median) %>% ungroup()

col_median <- intensity_row_median %>% group_by(msrun) %>% summarise(total_msrun_shift=median(norm_factor_row)) %>% ungroup()

pep_int_normalized <- msdata$pepmodstate_intensities %>%
  left_join(col_median) %>% mutate(norm_intensity=intensity*exp(-total_msrun_shift))

boxplotdata2 <-left_join(pep_int_normalized, select(msdata$msruns, msrun)) %>%
  mutate(log2intensity_norm = log2(norm_intensity),
         log2intensity = log2(intensity)) %>% 
  select(log2intensity, log2intensity_norm, msrun)

(p <- ggplot(boxplotdata2, aes(x=msrun, y=log2intensity)) +
    geom_boxplot() + theme(axis.text.x = element_text(angle = 90)))

msruns_hnorm <- list()
msruns_hnorm[["mschannel_shifts"]] <- msdata_full$msruns %>% left_join(col_median)
msdata$msrun_shifts <- msruns_hnorm[["mschannel_shifts"]]
# no need to estimate the effect scales per condition

rmsglmdata_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_data_', mstype, "_", data_info$fit_ver, '.RData'))
message('Saving data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata, msglm_def, msruns_hnorm,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', mstype, "_", data_info$data_ver, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full, fasta.dfs,
     file = rfulldata_filepath)

message('Done.')

