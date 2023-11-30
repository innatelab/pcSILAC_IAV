# assembly and analysis of Philipp's exp76 - GPKOW KD HeLa cells + SC35M (renilla)
# Experiments done by Philipp Hubel in 2015
# Author: Yiqi Huang
###############################################################################

project_id <- 'yhuang_iav'
message('Project ID=', project_id)
data_version <- "20231010"
fit_version <- "20231011b" 
mstype <- "phubel_exp123"
message('Dataset version is ', data_version)
message("Assembling fit results for project ", project_id,
        " (dataset v", data_version, ", fit v", fit_version, ")")

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))
require(msglm)

mop.max_nprocesses <- 32
mop.nprocesses <- 32
source(file.path(pipeline_scripts_path, 'init_cluster.R'))

require(dplyr)
require(stringr)
require(maxquantUtils)
require(furrr)

message('Loading data...')
load(file.path(scratch_path, paste0(project_id, '_msglm_data_', mstype,  "_", fit_version, '.RData')))
load(file.path(scratch_path, paste0(project_id, '_msdata_full_', mstype, "_", data_version, '.RData')))

message('Loading MSGLM model fit results...')

modelobj <- msdata$msentities[["object"]]
quantobj <- msdata$msentities[["quantobject"]]
modelobjs_df <- msdata$objects
modelobj_idcol <- paste0(modelobj, "_id")
quantobj_idcol <- paste0(quantobj, "_id")

fit_path <- file.path(scratch_path, paste0("msglm_", project_id, "_", mstype))#_', fit_version))
fit_files <- list.files(fit_path, paste0("msglm_", project_id, "_", mstype, "_", fit_version, '_\\d+\\.RData'))
message('Found ', length(fit_files), ' model file(s)')
fit_files.df <- tibble(filename = as.character(fit_files)) %>%
  tidyr::extract(filename, "chunk", ".+_(\\d+).RData$", convert=TRUE, remove=FALSE) %>%
  dplyr::left_join(dplyr::select(modelobjs_df, chunk, object_id), by="chunk") %>%
  dplyr::arrange(chunk)
#fit_diff <- anti_join(fit_files_new.df, fit_files.df)
require(RMySQL)
chunk_dispatcher_conn <- dbConnect(RMySQL::MySQL(),
                                   dbname="inlab_computing", user="inlab_dispatcher",
                                   host="tumevi4-websrv1.srv.mwn.de",
                                   password=Sys.getenv("INLAB_DISPATCHER_PASSWD"),
                                   timeout=300)
chunk_statuses.df <- dbGetQuery(chunk_dispatcher_conn,
                                str_c("SELECT * FROM jobchunks WHERE user='ge54heq2' AND job_id='",
                                      project_id, "_", mstype, "_", fit_version, "'")) %>%
  dplyr::mutate(start_time = as.POSIXlt(start_time),
                end_time = as.POSIXlt(end_time),
                fit_time = end_time - start_time)
dbDisconnect(chunk_dispatcher_conn)
table(chunk_statuses.df$status)
fit_files.df <- dplyr::inner_join(fit_files.df,
                                  dplyr::filter(chunk_statuses.df, status == "complete") %>%  #& fit_time > 0) %>%
                                    dplyr::select(chunk, status, fit_time), by="chunk")
# load fit results in parallel ----
plan(multisession, workers = 32)
fit_chunks <- seq_len(nrow(fit_files.df)) %>%
  furrr::future_map(.progress = TRUE, .options = furrr_options(stdout=FALSE, packages=c("msglm"), globals=c("fit_path", "fit_files.df")),
                    ~load_fit_chunk(.x))
names(fit_chunks) <- purrr::map_chr(fit_chunks, ~paste0(.$results_info$fit_version, '_',
                                                        .$msglm_results$objects$stats$object_id[1]))

fit_stats <- combine_fit_chunks(fit_chunks, 'stats')
fit_contrasts <- combine_fit_chunks(fit_chunks, 'contrast_stats')

rm(fit_chunks)

message('Done.')

if (modelobj == "protgroup") {
  # FIXME stats should be quantobj-dependent
  obj_conditions.df <- tidyr::expand(msdata_full$protgroup_intensities, protgroup_id, msrun) %>%
    dplyr::left_join(msdata_full$protgroup_intensities) %>%
    dplyr::inner_join(msdata$msruns) %>%
    dplyr::mutate(is_quanted = !is.na(intensity),
                  is_idented = replace_na(ident_type == "By MS/MS", FALSE)) %>%
    dplyr::group_by(condition, protgroup_id) %>%
    dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_quanted]),
                     nmsruns_idented = n_distinct(msrun[is_idented])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(object_id = protgroup_id)
  
  modelobjs_df <- dplyr::mutate(modelobjs_df,
                                is_msvalid_object = npepmods_unique >= 2L)#(nprotgroups_sharing_proteins == 1 || nproteins_have_razor > 0))
} else if (modelobj == "protregroup") {
  obj_conditions.df <- expand(msdata_full$pepmodstate_intensities, msrun, pepmod_id) %>%
    dplyr::inner_join(msdata$msruns) %>%
    dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
    dplyr::mutate(is_quanted = !is.na(intensity),
                  is_idented = is_quanted) %>%
    dplyr::left_join(filter(msdata$protregroup2pepmod, is_specific)) %>%
    dplyr::group_by(condition, protregroup_id) %>%
    dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_idented]), # count msruns
                     nmsruns_idented = n_distinct(msrun[is_quanted])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(object_id = protregroup_id)
  
  modelobjs_df <- dplyr::mutate(modelobjs_df,
                                is_msvalid_object = npepmods_unique >= 2L)
}

contrastXcondition.df <- as_tibble(as.table(msglm_def$conditionXmetacondition)) %>% dplyr::filter(n != 0) %>%
  dplyr::select(-n) %>%
  dplyr::inner_join(as_tibble(as.table(msglm_def$metaconditionXcontrast))) %>% dplyr::filter(n != 0) %>% 
  dplyr::arrange(contrast, metacondition, condition)

pre_object_contrasts.df <- obj_conditions.df %>% 
  group_by(object_id, condition) %>%
  mutate(obs = seq(n())) %>%
  ungroup %>%
  # complete the data
  tidyr::complete(object_id, condition, obs) %>%
  select(-obs) %>%
  mutate(protregroup_id = coalesce(protregroup_id, object_id)) %>% 
  replace_na(list( nmsruns_quanted = 0, nmsruns_idented = 0)) %>% 
  dplyr::inner_join(contrastXcondition.df) %>%
  dplyr::mutate(is_lhs = n > 0) %>%
  dplyr::group_by(object_id, contrast, is_lhs) %>%
  dplyr::summarise(has_quanted = any(!is.na(nmsruns_quanted)),
                   nmsruns_quanted_min = min(nmsruns_quanted, na.rm=TRUE),
                   nmsruns_quanted_max = max(nmsruns_quanted, na.rm=TRUE),
                   has_idented = any(!is.na(nmsruns_idented)),
                   nmsruns_idented_min = min(nmsruns_idented, na.rm=TRUE),
                   nmsruns_idented_max = max(nmsruns_idented, na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(nmsruns_quanted_min = if_else(has_quanted, nmsruns_quanted_min, 0L),
                nmsruns_quanted_max = if_else(has_quanted, nmsruns_quanted_max, 0L),
                nmsruns_idented_min = if_else(has_idented, nmsruns_idented_min, 0L),
                nmsruns_idented_max = if_else(has_idented, nmsruns_idented_max, 0L)) %>%
  dplyr::group_by(object_id, contrast) %>%
  dplyr::summarise(nmsruns_quanted_lhs_min = nmsruns_quanted_min[is_lhs],
                   nmsruns_quanted_lhs_max = nmsruns_quanted_max[is_lhs],
                   nmsruns_idented_lhs_min = nmsruns_idented_min[is_lhs],
                   nmsruns_idented_lhs_max = nmsruns_idented_max[is_lhs],
                   nmsruns_quanted_rhs_min = nmsruns_quanted_min[!is_lhs],
                   nmsruns_quanted_rhs_max = nmsruns_quanted_max[!is_lhs],
                   nmsruns_idented_rhs_min = nmsruns_idented_min[!is_lhs],
                   nmsruns_idented_rhs_max = nmsruns_idented_max[!is_lhs]) %>%
  dplyr::ungroup() # In our case, the condition and metacondition is the same, so the nmsruns_quanted/idented will always have the same min and max, but it's not always the case if there're multiple conditions in one metacondition (e.g. in AP-MS analysis)

contrasts.df <- dplyr::ungroup(msglm_def$contrasts) %>%
  dplyr::inner_join(tidyr::pivot_wider(dplyr::mutate(as.data.frame.table(msglm_def$metaconditionXcontrast, responseName="w"),
                                                     side = if_else(w > 0, "lhs", "rhs")) %>% dplyr::filter(w != 0),
                                       c(contrast), names_from = "side", values_from = "metacondition",
                                       names_glue = "{.value}_{side}")) %>%
  dplyr::mutate(offset = 0, offset_prior = 0) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_lhs = condition, ko_lhs = ko, treatment_lhs = treatment)) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_rhs = condition, ko_rhs = ko, treatment_rhs = treatment))

object_contrasts_thresholds.df <- contrasts.df %>%
  mutate(p_value_threshold = case_when(TRUE ~ 0.001),
         #p_value_threshold_lesser = case_when(TRUE ~ 0.01),
         median_threshold = case_when(TRUE ~ 0.5)#,
         #median_threshold_lesser = case_when(TRUE ~ 0.125)
  )

object_contrasts.df <- fit_contrasts$object_conditions %>%
  inner_join(pre_object_contrasts.df) %>% 
  dplyr::filter(var %in% c("obj_cond_labu", "obj_cond_labu_replCI")) %>%
  dplyr::inner_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_valid_comparison = (pmax(nmsruns_quanted_lhs_max, nmsruns_quanted_rhs_max)>=3), #quantified in 3/4 on either side
                is_signif = p_value <= p_value_threshold & abs(median) >= median_threshold,
                #is_signif_lesser = p_value <= p_value_threshold_lesser & abs(median) >= median_threshold_lesser,
                is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse,
                is_hit = is_hit_nomschecks & is_valid_comparison, 
                change = if_else(is_signif, if_else(median < 0, "-", "+"), ".")) 

object_contrast_stats.df <- dplyr::group_by(object_contrasts.df, contrast, contrast_type, ci_target) %>%
  dplyr::summarise(p_value_001 = quantile(p_value, 0.001),
                   p_value_01 = quantile(p_value, 0.01),
                   p_value_05 = quantile(p_value, 0.05),
                   median_abs_50 = quantile(abs(median[p_value <= 0.1]), 0.5),
                   median_abs_95 = quantile(abs(median[p_value <= 0.1]), 0.95),
                   median_abs_99 = quantile(abs(median[p_value <= 0.1]), 0.99),
                   n_hits = sum(is_hit_nomschecks, na.rm = TRUE),
                   n_plus = sum(change == "+"),
                   n_minus = sum(change == "-")) %>%
  dplyr::ungroup() %>% 
  filter(ci_target == "average")

View(filter(object_contrast_stats.df, ci_target == "average") %>% dplyr::arrange(desc(n_hits)))

object_contrasts_wide.df <- pivot_wider(object_contrasts.df,
                                        id_cols = c("ci_target", "object_id", "object_label", "is_viral"),
                                        names_from = "contrast", names_sep = ".",
                                        values_from = c("median", "mean", "sd", "p_value", "is_hit", "change"))

rfit_filepath <- file.path(results_path, paste0(project_id, '_msglm_fit_', mstype, '_', fit_version, '.RData'))
results_info <- list(project_id = project_id, mstype=mstype,
                     data_version = data_version, fit_version = fit_version,
                     modelobj = modelobj, quantobj = quantobj)
message('Saving full analysis results to ', rfit_filepath, '...')
save(results_info, fit_stats, fit_contrasts,
     object_contrasts.df, object_contrasts_wide.df,
     object_contrasts_thresholds.df,
     file = rfit_filepath)
message('Done.')

report_cols <- c("object_label", "object_id", "gene_names",
                 "majority_protein_acs", "protein_descriptions",
                 "is_contaminant", "is_viral")

objects4report.df <- dplyr::select(msdata$objects, any_of(report_cols)) %>%
  dplyr::semi_join(dplyr::select(object_contrasts.df, object_id))

object_contrasts_report.df <- objects4report.df %>%
  dplyr::left_join(pivot_wider(object_contrasts.df, c(ci_target, object_id),
                               names_from = "contrast", values_from = c("change", "is_valid_comparison", "is_hit",  "median", "p_value", "sd"),
                               names_sep=".")) %>%
  dplyr::arrange(gene_names, majority_protein_acs, ci_target)

write_tsv(object_contrasts_report.df,
          file.path(analysis_path, "reports", paste0(project_id, '_', mstype, '_contrasts_report_', fit_version, '_wide.txt')))

full_gene_list <- object_contrasts_report.df %>% 
  filter(ci_target == "average") %>% 
  select(gene_names) %>% 
  mutate(gene_names = str_remove_all(gene_names, ";.*"))
write_tsv(full_gene_list,
          file.path(analysis_path, "reports", paste0(project_id, '_', mstype, '_all_genes_', fit_version, '.txt')))

#make supplementary table for publication ----
load(rfit_filepath)

