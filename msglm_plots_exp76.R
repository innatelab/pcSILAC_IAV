# plotting the volcanos for contrasts
# peptide heatmaps
# Yiqi Huang

#The following part is only needed if starting from a fresh environment----
project_id <- 'yhuang_iav'
message('Project ID=', project_id)
data_version <- "20231010"
fit_version <- "20231011" 
mstype <- "phubel_exp76"
msfolder <- "phubel_exp76+123"
message("Project ID=", project_id, " data version=", data_version)

require(tidyverse)
require(rlang)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

message('Loading data...') 
load(file.path(scratch_path, str_c(project_id, '_msdata_full_', mstype, "_", data_version, '.RData')))
load(file.path(scratch_path, str_c(project_id, '_msglm_data_', mstype,  "_", fit_version, '.RData')))
#load(file.path(scratch_path, str_c(project_id, '_msglm_fit_meanfield_', mstype, "_", fit_version, '.RData')))
load(file.path(results_path, str_c(project_id, '_msglm_fit_', mstype, "_", fit_version, '.RData')))

modelobj <- msdata$msentities['object']
quantobj <- msdata$msentities['quantobject']
modelobjs_df <- msdata$objects
modelobj_idcol <- paste0(modelobj, "_id")
quantobj_idcol <- paste0(quantobj, "_id")

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
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_lhs = condition, kd_lhs = kd)) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_rhs = condition, kd_rhs = kd))

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

#plotting starts from here----
require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)
require(ggpubr)

obj_labu_shift <- msdata[[str_c(quantobj, "_mscalib")]]$zShift
modelobj_suffix <- "_exp76"

source(file.path(misc_scripts_path, 'ggplot_ext.R'))
source(file.path(misc_scripts_path, 'furrr_utils.R'))

kd_palette <- c("Scrambled" = "gray", "siGPKOW" = "#F9CB40", "siJAK1" = "red")
hit_palette <- c("non-hit" = "grey", "hit" = "black", "viral hit" = "#F9CB40", "only sig" = "light blue")
base_font_family <- "Segoe UI Symbol"
base_plot_path <- file.path(analysis_path, 'plots', str_c(mstype, "_",data_version, "_", fit_version))
sel_ci_target <- "average"


object_contrasts_thresholds.df <- dplyr::select(contrasts.df, offset, offset_prior, contrast, contrast_type) %>%
  dplyr::mutate(
    p_value_threshold = case_when(contrast_type == "filtering" ~ 1E-3,
                                  contrast_type == "comparison" ~ 1E-3,
                                  TRUE ~ NA_real_),
    median_threshold = case_when(contrast_type == "filtering" ~ pmax(2.0, 2.0 + abs(offset - offset_prior)),
                                 contrast_type == "comparison" ~ pmax(0.5, 0.25 + abs(offset - offset_prior)),
                                 TRUE ~ NA_real_),
    median_max = case_when(contrast_type == "filtering" ~ 12,
                           contrast_type == "comparison" ~ 6,
                           TRUE ~ NA_real_)
  )

#volcano for all contrasts----
object_contrasts_4show.df <- object_contrasts.df %>%
  dplyr::filter(var %in% c("obj_cond_labu", "obj_cond_labu_replCI")) %>%
  dplyr::ungroup() %>%
  select(-contains("threshold")) %>%
  dplyr::inner_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - offset) >= median_threshold),
                is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse,
                is_hit = is_hit_nomschecks & (pmax(nmsruns_quanted_lhs_max, nmsruns_quanted_rhs_max)>=3), #quantified in 3/4 on either side
                hit_type = case_when(is_hit & is_viral ~ "viral hit",
                                     is_hit ~ "hit", 
                                     is_hit_nomschecks ~ "only sig",
                                     TRUE ~ "non-hit"), 
                is_sel = is_hit_nomschecks,
                mean_trunc = pmax(-median_max, pmin(median_max, mean - offset)) + offset,
                median_trunc = pmax(-median_max, pmin(median_max, median - offset)) + offset,
                truncation = scatter_truncation(median, median_trunc, -log10(p_value), -log10(p_value),
                                                is_hit_nomschecks | !is_signif),
                truncation_type = point_truncation_type(truncation, is_signif),
                show_label = coalesce(is_hit_nomschecks, FALSE),
                short_label = str_remove_all(object_label, "Oxidation_|Phospho_|_M\\d$")) %>%
  dplyr::group_by(contrast, ci_target) %>%
  tidyr::separate(contrast, into = c("group1", "group2"), sep = "_vs_", remove = FALSE) %>%
  dplyr::mutate(show_label = if_else(rep.int(sum(show_label & median > 0) >= 400L, n()), is_hit, show_label)) %>% #what does this do?
  dplyr::ungroup()

require(furrr)
plan(multicore, workers = 32)
plot_furrr_opts <- furrr_options(globals = c("base_plot_path", "base_font_family",
                                             "project_id", "data_version", "fit_version", "obj_labu_shift",
                                             "mlog10_trans", "mlog_pow_trans", "mlog_breaks",
                                             "theme_bw_ast",
                                             "point_truncation_shape_palette", "point_truncation_size_palette",
                                             "treatment_palette",
                                             "sel_ci_target", "modelobj", "quantobj", "modelobj_idcol", "quantobj_idcol", "modelobjs_df",
                                             "msglm_def", "fit_stats", "msdata", "msdata_full", "object_contrasts_4show.df"),
                                 packages = c("dplyr", "ggplot2", "Cairo", "ggrepel", "ggrastr", "stringr"),
                                 stdout = TRUE)

group_by(object_contrasts_4show.df, ci_target, contrast,
         offset, median_threshold, p_value_threshold) %>%
  group_walk(.keep = TRUE,
             #future_group_walk(.progress=TRUE, .keep=TRUE, .options=plot_furrr_opts,
             function(sel_object_contrast.df, contrast_info) {
               message("Plotting ", contrast_info$contrast, " ci_target=", contrast_info$ci_target)
               labels_lhs.df <- dplyr::filter(sel_object_contrast.df, (median > offset) & (is_sel | is_signif & show_label))
               if (nrow(labels_lhs.df) > 300) {
                 labels_lhs.df <- dplyr::filter(labels_lhs.df, is_sel)
               }
               labels_rhs.df <- dplyr::filter(sel_object_contrast.df, (median < offset) & (is_sel | is_signif & show_label))
               if (nrow(labels_rhs.df) > 300) {
                 labels_rhs.df <- dplyr::filter(labels_rhs.df, is_sel)
               }
               labels.df <- bind_rows(labels_lhs.df, labels_rhs.df)
               nlabels <- nrow(labels.df)
               
               p <- ggplot(sel_object_contrast.df,
                           aes(x=median_trunc, y=p_value, shape=truncation, size=truncation_type, color=hit_type)) +
                 geom_hline(data=contrast_info, aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
                 #geom_hline(data=contrast_info, aes(yintercept = p_value_max), linetype=1, color="darkgray") +
                 geom_vline(data=contrast_info, aes(xintercept = offset), linetype=1, color="darkgray") +
                 geom_vline(data=contrast_info, aes(xintercept = offset + median_threshold), linetype=2, color="darkgray") +
                 geom_vline(data=contrast_info, aes(xintercept = offset - median_threshold), linetype=2, color="darkgray") +
                 geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_signif),
                                 alpha = 0.1, size = 0.5, color = "darkgray") +
                 geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & !is_hit), shape=1) +
                 geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & is_hit)) +
                 geom_text_repel(data = labels.df,
                                 aes(label = short_label),
                                 size = if_else(nlabels > 20, 2.5, 3.5),
                                 force = if_else(nlabels > 20, 0.25, 1.0),
                                 nudge_y = -0.12,
                                 point.padding = 0.02,
                                 box.padding = if_else(nlabels > 20, 0.1, 0.25),
                                 show.legend = FALSE, segment.color = "gray",
                                 max.overlaps = Inf) +
                 scale_y_continuous(trans = mlog10_trans(), limits = c(1.0, NA)) +
                 #scale_fill_gradient(low="gray75", high="black") +
                 #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
                 scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
                 scale_size_manual(values=point_truncation_size_palette, guide="none") +
                 scale_color_manual(values = hit_palette, na.value = "magenta") +
                 #facet_grid(p_value_range ~ contrast, scales = "free_y") +
                 ggtitle(contrast_info$contrast, subtitle=str_c("ci_target=", contrast_info$ci_target)) +
                 theme_bw_ast(base_family = base_font_family)
               plot_path <- file.path(base_plot_path, "meanfield", str_c("volcanos_contrasts_", contrast_info$ci_target))
               if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
               
               ggsave(filename = file.path(plot_path, 
                                           str_c(project_id, '_', fit_version, '_volcano_',
                                                 str_replace_all(contrast_info$contrast, ":|@", "_"), '.pdf')),
                      plot=p, width=15, height=18, device=cairo_pdf, family=base_font_family) #For windows users only: beware of what font you're using! If it's the first time you use this code, you need to import the fonts with "extrafont" package
             })

#Making timecourse for all proteins. The dots are from the original LFQ intensity values of MQ.----
#sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(gene_names, "IFIT1"))  #This is used for debugging
#sel_objects.df <- dplyr::semi_join(modelobjs_df, dplyr::select(fit_diff, object_id), by = "object_id")
sel_objects.df <- dplyr::semi_join(modelobjs_df, dplyr::select(fit_stats$objects, object_id), by="object_id")#This is all!

dplyr::group_by(sel_objects.df, object_id) %>% do({
  sel_obj.df <- .
  message("Plotting ", sel_obj.df$object_label, " time course")
  
  sel_var <- if_else(sel_ci_target == "average", "obj_cond_labu", "obj_cond_labu_replCI")
  sel_obj_conds.df <- dplyr::inner_join(dplyr::select(sel_obj.df, object_id),
                                        dplyr::filter(fit_stats$object_conditions, var == sel_var)) %>%
    dplyr::inner_join(dplyr::select(msglm_def$conditions, condition, treatment, timepoint, timepoint_num), by="condition") %>%
    dplyr::arrange(timepoint_num) %>% 
    dplyr::mutate(dplyr::across(c(mean, median, starts_with("q")),
                                ~2^(.x + obj_labu_shift))) %>%
    dplyr::mutate(q97.5_limit = max(median + (q75 - q25)*5))
  
  sel_obj_contrasts.df <- dplyr::semi_join(dplyr::filter(object_contrasts_4show.df, ci_target == sel_ci_target & str_starts(var, "obj_cond_labu"),
                                                         group1 == "MPXV"),
                                           sel_obj.df, by="object_id") %>% 
    mutate(y.position = max(sel_obj_conds.df$`q97.5`)*1.1,
           treatment = group1,
           timepoint = as.numeric(str_extract(group2, "(?<=@)\\d+")),
           #p_value = formatC(p_value, format = "e", digits = 2),
           p_label = case_when(p_value < 0.001 ~ "***",
                               p_value < 0.01 ~ "**",
                               p_value < 0.05 ~ "*",
                               TRUE ~ ""))# for the pvalue labels in the graph
  if (modelobj == "protgroup") {
    sel_obj_msdata.df <- sel_obj.df %>%
      dplyr::inner_join(msdata$protgroup_intensities) %>%
      dplyr::inner_join(msdata$msruns)
  } else if (modelobj == "protregroup") {
    sel_obj_msdata.df <- sel_obj.df %>%
      dplyr::inner_join(msdata_full$protein2protregroup, by = "protregroup_id") %>%
      dplyr::group_by(object_id) %>%
      dplyr::filter(is_majority == any(is_majority)) %>%
      dplyr::ungroup() %>%
      dplyr::select(object_id, protregroup_id, protein_ac) %>%
      dplyr::inner_join(msdata_full$protein2protgroup) %>%
      dplyr::group_by(object_id, protgroup_id) %>%
      dplyr::filter(is_majority == any(is_majority)) %>%
      dplyr::ungroup() %>%
      dplyr::select(object_id, protgroup_id, protregroup_id) %>% dplyr::distinct() %>%
      dplyr::inner_join(msdata_full$protgroup_intensities) %>%
      dplyr::right_join(msdata$msruns, by="msexperiment_mq") %>%
      dplyr::left_join(dplyr::select(msdata$msrun_shifts, msrun, total_msrun_shift), by="msrun") %>%
      dplyr::mutate(intensity_norm = intensity*2^(-total_msrun_shift)) %>%
      dplyr::inner_join(msdata_full$msruns)
    #print(sel_obj_msdata.df)
    obj_shifts.df = dplyr::inner_join(sel_obj_msdata.df,
                                      dplyr::select(sel_obj_conds.df, condition, median), by="condition") %>%
      dplyr::group_by(object_id) %>%
      dplyr::summarise(obj_shift = median(log2(intensity_norm) - log2(median), na.rm=TRUE), .groups = "drop")
    #print(obj_shifts.df)
    sel_obj_msdata.df <- dplyr::left_join(sel_obj_msdata.df, obj_shifts.df) %>%
      dplyr::mutate(intensity_norm = intensity_norm * 2^(-obj_shift))
  }
  
  #print(sel_obj_msdata.df)
  if (nrow(sel_obj_conds.df) > 0) {
    p <-
      ggplot(data=sel_obj_conds.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
      geom_ribbon(aes(x = timepoint_num, ymin = `q2.5`, ymax=`q97.5`),
                  alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
      geom_ribbon(aes(x = timepoint_num, ymin = q25, ymax= q75),
                  alpha=0.5, stat = "identity", size=0.5) +
      geom_path(aes(x = timepoint_num, y = median), alpha=0.5, size=1, stat="identity") +
      geom_point(data=sel_obj_msdata.df,
                 aes(y = intensity_norm),
                 position = position_jitter(width = 0.75, height = 0, seed=12323), size=1) +
      geom_text(data = sel_obj_contrasts.df, aes(x = timepoint, y = y.position, label = p_label), colour = "black")+
      theme_bw_ast(base_family = "", base_size = 8) +
      scale_x_continuous(breaks=unique(msdata$msruns$timepoint_num)) +
      scale_color_manual(values=treatment_palette) +
      scale_fill_manual(values=treatment_palette) +
      scale_y_log10() +
      ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
              subtitle=sel_obj.df$protein_description) +
      facet_wrap( ~ object_label, ncol =1, scales = "free")
    plot_path <- file.path(base_plot_path,
                           str_c("timecourse_", sel_ci_target,
                                 modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")))
    if (!dir.exists(plot_path)) {dir.create(plot_path, recursive = TRUE)}
    ggsave(p, file = file.path(plot_path, str_c(project_id, "_", data_info$msfolder, '_', fit_version, "_",
                                                str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], ".pdf")),
           width=8, height=6, device = cairo_pdf)
  }
  tibble()
})

#Making peptide heatmaps for every protein group----
#sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(gene_names, "IFIT1"))
sel_objects.df <- modelobjs_df # all!!!

sel_pepmodstates.df <- dplyr::inner_join(sel_objects.df,
                                         msdata_full[[str_c(modelobj, "2pepmod")]]) %>% 
  dplyr::inner_join(msdata_full$pepmodstates) %>%
  dplyr::inner_join(dplyr::transmute(msdata_full$pepmodstate_intensities, pepmodstate_id,
                                     is_idented,
                                     is_quanted = !is.na(intensity)) %>%
                      dplyr::group_by(pepmodstate_id) %>%
                      dplyr::summarise(nidented = sum(is_idented, na.rm = TRUE),
                                       nquanted = sum(is_quanted, na.rm = TRUE),
                                       .groups = "drop")) %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmods, pepmod_id, peptide_id, pepmod_seq, peptide_seq)) %>%
  dplyr::select(object_id, object_label, pepmod_id, majority_protein_acs, protac_label, gene_label, gene_names, #msfraction,
                nidented, nquanted, pepmod_seq, peptide_seq, charge, pepmodstate_id, is_specific)

if (exists("fit_stats") && has_name(fit_stats, "quantobjects")) {
  sel_pepmodstates.df <- dplyr::left_join(sel_pepmodstates.df,
                                          dplyr::select(dplyr::filter(fit_stats$quantobjects, var=="qobj_shift"),
                                                        pepmodstate_id=quantobject_id, pms_median = median),
                                          by="pepmodstate_id") %>%
    dplyr::arrange(object_id, desc(is_specific), desc(tidyr::replace_na(pms_median, -1000.0)), desc(nidented), desc(nquanted),
                   peptide_seq, pepmod_seq, charge)
} else {
  sel_pepmodstates.df <- dplyr::arrange(sel_pepmodstates.df, object_id, desc(is_specific), desc(nidented), desc(nquanted),
                                        peptide_seq, pepmod_seq, msfraction, charge) %>%
    mutate(pms_median = NA_real_)
}
sel_pepmodstates.df <- dplyr::mutate(sel_pepmodstates.df,
                                     # FIXME remove pepmod_id once mod_seq is really mod_seq
                                     pepmodstate_ext = paste0(pepmod_seq, ".", charge, #" F", msfraction,
                                                              " (", pepmodstate_id, ") ",
                                                              if_else(is.na(pms_median),
                                                                      if_else(is_specific, "+", ""), "*"))
                                     %>% factor(., levels = unique(.)))

sel_pepmod_intens.df <- tidyr::crossing(pepmodstate_id = sel_pepmodstates.df$pepmodstate_id,
                                        msrun = unique(msdata_full$pepmodstate_intensities$msrun)) %>%
  dplyr::inner_join(dplyr::select(sel_pepmodstates.df, object_id, object_label, 
                                  pepmod_id, pepmodstate_id, charge, #msfraction, 
                                  pepmodstate_ext, pms_median) %>%
                      dplyr::distinct()) %>%
  dplyr::inner_join(distinct(select(msdata$msruns, msexperiment, msfraction, msrun, condition, kd, replicate)) %>%
                      dplyr::left_join(dplyr::select(msdata$msrun_shifts, msrun, total_msrun_shift), by="msrun") %>%
                      dplyr::arrange(kd, replicate, 
                                     msrun) %>%
                      dplyr::mutate(msrun.2 = factor(msrun, levels=unique(msrun)),
                                    msexperiment.2 = factor(msexperiment, levels=unique(msexperiment))))%>%
  dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
  dplyr::mutate(intensity_norm = 2^(-total_msrun_shift)*intensity,
                msrun = msrun.2, msrun.2 = NULL,
                msexperiment = msexperiment.2, msexperiment.2 = NULL) %>%
  dplyr::group_by(pepmodstate_id) %>%
  dplyr::mutate(intensity_norm_trunc = pmin(intensity_norm, quantile(intensity_norm, 0.95, na.rm=TRUE))) %>%
  dplyr::ungroup()

sel_pepmods.df <- dplyr::group_by(sel_pepmod_intens.df, pepmod_id) %>%
  dplyr::summarize(n_pepmod_quants = sum(!is.na(intensity)),
                   median_quant = median(intensity, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_pepmod_quants), desc(median_quant), pepmod_id)

group_by(sel_pepmod_intens.df, object_id) %>%
  #group_walk(.keep=TRUE,
  future_group_walk(.progress = TRUE, .keep = TRUE, .options = plot_furrr_opts,
                    function(obj_pepmod_intens.df, obj_row) {
                      obj_row <- dplyr::semi_join(modelobjs_df, obj_row, by="object_id")
                      shown_pepmod_intens.df <- mutate(obj_pepmod_intens.df, !!modelobj_idcol := object_id,
                                                       intensity_norm = pmax(pmin(intensity_norm, quantile(intensity_norm, 0.95, na.rm=TRUE)),
                                                                             quantile(intensity_norm, 0.05, na.rm=TRUE))) #sel_pepmod_intens.df
                      obj_label <- str_remove(obj_row$object_label, "\\.\\.\\.$")
                      message("Plotting ", obj_label, "...")
                      p <- ggplot(shown_pepmod_intens.df) +
                        geom_tile(aes(x = msexperiment, y = pepmodstate_ext,
                                      fill = intensity_norm, color = ident_type),
                                  na.rm = FALSE, size=0.5, width=0.85, height=0.85) +
                        theme_bw_ast(base_family = base_font_family, base_size = 10) +
                        theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0)) +
                        guides(color=guide_legend("ident_type", override.aes = list(fill=NA, size=2))) +
                        ggtitle(str_c(obj_label,  " (", modelobj_idcol, "=", obj_row$object_id,
                                      ", ac=", obj_row$protac_label, ") peptide map"),
                                subtitle = obj_row$protein_description) +
                        scale_fill_distiller(na.value="#00000000", type="div", palette = "Spectral") +
                        scale_y_discrete(breaks = levels(obj_pepmod_intens.df$pepmodstate_ext)) +
                        scale_color_manual(na.value="#00000000",
                                           values=c("MULTI-MSMS"="black", "MULTI-MATCH-MSMS"="khaki",
                                                    "ISO-MSMS"="slateblue",
                                                    "MSMS"="cornflowerblue", "MULTI-SECPEP"="firebrick",
                                                    "MULTI-MATCH"="gray"))
                      
                      plot_path <- file.path(base_plot_path, "pepmodstate_heatmaps")
                      if (!dir.exists(plot_path)) dir.create(plot_path)
                      ggsave(p, file = file.path(plot_path,
                                                 paste0(project_id, "_", data_version, "_pepmodstate_heatmap_",
                                                        obj_label, "_", obj_row$object_id, ".pdf")),
                             width=2 + 0.5*n_distinct(shown_pepmod_intens.df$msrun),
                             height=2 + 2*n_distinct(shown_pepmod_intens.df$object_id) +
                               1.5*min(20, 0.12*n_distinct(shown_pepmod_intens.df$pepmodstate_id)),
                             device=cairo_pdf, family=base_font_family, limitsize=FALSE)
                      gc()
                    })