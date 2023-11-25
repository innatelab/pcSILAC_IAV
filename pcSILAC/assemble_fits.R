# assembly and analysis of pulse dynamics fits
# 
# Author: Alexey Stukalov
###############################################################################

print(commandArgs())

job.args <- c("phubel_pulse", "20171025", "0")
if (!exists('job.args')) {
    job.args <- commandArgs(trailingOnly = TRUE)
}

project_id <- job.args[[1]]
message('Project ID=', project_id)

job_version <- job.args[[2]]
message('Job dataset version is ', job_version)

job_id <- as.integer(job.args[[3]])
message('Job ID is ', job_id)

source('~/R/config.R')
source(file.path(base_scripts_path, 'misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))
source(file.path(misc_scripts_path, 'reshape.R'))
source(file.path(base_scripts_path, 'misc/ggplot_ext.R'))
source(file.path(base_scripts_path, "msglm/R/frame_utils.R"))
source(file.path(base_scripts_path, 'msglm/R/assemble_utils.R') )

analysis_path <- scratch_path

mop.max_nprocesses <- 8
mop.nprocesses <- 8
source(file.path(pipeline_scripts_path, 'init_cluster.R'))

require(CeMMmisc)
require(rstan)
require(dplyr)
require(rjson)

message('Loading data...')
load(file.path(scratch_path, paste0(project_id, "_", job_version, "_lean.RData")))

strip_samples <- TRUE

message( 'Loading pulse dynamics models...' )

analysis_path <- scratch_path
fit_path <- file.path(analysis_path, paste0("dynamics_fit_rescaled_", job_version))
fit_files <- list.files(fit_path, paste0('pulse_dynamics_fit_', job_version, '_\\d+\\.RData'))
message( 'Found ', length( fit_files ), ' model file(s)' )
fit_files.df <- bind_rows(lapply(strsplit(as.character(fit_files), '[_.]', fixed=FALSE),
                                        function( file_chunks ) {
  data.frame(protgroup_id = as.integer(file_chunks[[5]]),
             stringsAsFactors = FALSE) } )) %>%
  dplyr::mutate(filename = as.character(fit_files),
                task_id = protgroup_id + 1,
                date = file.mtime(file.path(fit_path, filename))) %>%
  dplyr::arrange(task_id)
id_range_breaks <- which(c(fit_files.df$task_id[-1] - 1L, sum(dplyr::mutate(ms_data$protgroups, is_fit=TRUE)$is_fit)) != 
                           c(fit_files.df$task_id[-length(fit_files.df$task_id)], sum(dplyr::mutate(ms_data$protgroups, is_fit=TRUE)$is_fit)))
View(fit_files.df[sort(c(id_range_breaks, id_range_breaks+1L)), ] %>%
     dplyr::inner_join(dplyr::select(ms_data$protgroups, protgroup_id, majority_protein_acs, gene_names)))
View(dplyr::anti_join(dplyr::select(ms_data$protgroups, protgroup_id, majority_protein_acs, gene_names), fit_files.df))
write(paste0(setdiff(ms_data$protgroups$protgroup_id+1L, fit_files.df$task_id), collapse="\n"),
      file=file.path("~/projects/adhoc/phubel_pulse", paste0("missing_task_ids_", job_version)))

#load(file.path(fit_path, fit_files[1]))

if ( !is.null( mop.cluster ) ) {
  clusterEvalQ(mop.cluster, library(dplyr))
  clusterExport(mop.cluster, varlist=c('fit_path', 'fit_files.df', 'process_msglm_chunk'))
  fit_reports <- clusterApplyLB(mop.cluster, seq_along(fit_files.df$task_id), process_msglm_chunk)
  clusterEvalQ( mop.cluster, gc() )
} else {
  fit_reports <- lapply(seq_along(fit_files.df$task_id)[1:10], process_msglm_chunk)
}
names(fit_reports) <- sapply(fit_reports, function(report) paste0("pulse_fit_", report$analysis_info$protgroup_id))
# add protgroup information to each statistics frame
fit_reports <- lapply(fit_reports, function(report) {
  pg_id <- report$analysis_info$protgroup_id
  report$vars_results <- lapply(report$vars_results, function(res) {
    for (frname in c("stats", "condition_contrast_stats", "degradation_contrast_stats")) {
      if (frname %in% names(res)) {
        res[[frname]] <- dplyr::mutate(res[[frname]], protgroup_id = pg_id)
      }
    }
    return(res)
  })
  return(report)
})

fit_stats <- join_msglm_reports_allsections(fit_reports, "stats", results_tag="vars_results")

fit_cond_contrast_stats <- join_msglm_reports_allsections(fit_reports, "condition_contrast_stats", results_tag="vars_results")
fit_cond_contrast_stats <- fit_cond_contrast_stats[!sapply(fit_cond_contrast_stats, function(stats) is.null(stats) || nrow(stats)==0)]

fit_degr_contrast_stats <- join_msglm_reports_allsections(fit_reports, "degradation_contrast_stats", results_tag="vars_results")
fit_degr_contrast_stats <- fit_degr_contrast_stats[!sapply(fit_degr_contrast_stats, function(stats) is.null(stats) || nrow(stats)==0)]

fit_status.df <- bind_rows(lapply(fit_reports, function(report){
  data.frame(mcmc_converged = report$fit_status$mcmc_converged,
             mcmc_start_time = report$fit_status$mcmc_start_time,
             mcmc_end_time = report$fit_status$mcmc_end_time,
             job_chunk = report$analysis_info$job_chunk,
             protgroup_id = report$analysis_info$protgroup_id,
             gene_name = report$analysis_info$gene_name,
             stringsAsFactors = FALSE)
})) %>% dplyr::inner_join(dplyr::select(ms_data$protgroups, protgroup_id, is_full_quant, is_contaminant, is_reverse))

ms_data$protgroups <- dplyr::inner_join(ms_data$protgroups, dplyr::select(ms_data$protgroup_stats, protgroup_id, n_quanted)) %>%
  dplyr::arrange(protgroup_id)

exclude_protgroup_ids <- dplyr::filter(ms_data$protgroups, is_contaminant | is_reverse | (n_quanted == 0))$protgroup_id

# adjust P-values (might be wrong if there are multiple versions of the analysis for the same data)
for (df_name in names(fit_stats)) {
  df <- fit_stats[[df_name]]
  if ('p_value' %in% colnames(df)) {
    message('Adjusting P-values for ', df_name, '...')
    group_vars <- intersect(c('var', 'param', 'param_type'), colnames(df))
    fit_stats[[df_name]] <- df %>%
      group_by_(.dots = group_vars) %>%
      mutate(p_value_adj = p.adjust(p_value, method = "BY")) %>%
      ungroup() %>%
      mutate(is_hit = p_value_adj <= 1E-2 & !(protgroup_id %in% exclude_protgroup_ids))
  }
}
for (df_name in names(fit_cond_contrast_stats)) {
  message('Adjusting P-values for ', df_name, '...')
  df <- fit_cond_contrast_stats[[df_name]]
  group_vars <- intersect(c('var', 'param', 'index_contrast', 'contrast'),
                          colnames(df))
  df <- group_by_(df, .dots = group_vars) %>%
    mutate(p_value = 2*pmin(prob_nonneg, prob_nonpos),
           p_value_adj = p.adjust(p_value, method = "BY"),
           prob_nonneg_adj = p.adjust(prob_nonneg, method = "BY"),
           prob_nonpos_adj = p.adjust(prob_nonpos, method = "BY")) %>%
    ungroup() %>%
    mutate(is_hit = p_value_adj <= 1E-2 & !(protgroup_id %in% exclude_protgroup_ids))
  fit_cond_contrast_stats[[df_name]] <- df
}
for (df_name in names(fit_degr_contrast_stats)) {
  message('Adjusting P-values for ', df_name, '...')
  df <- fit_degr_contrast_stats[[df_name]]
  group_vars <- intersect(c('var', 'param', 'index_contrast', 'contrast'),
                          colnames(df))
  df <- group_by_(df, .dots = group_vars) %>%
    mutate(p_value = 2*pmin(prob_nonneg, prob_nonpos),
           p_value_adj = p.adjust( p_value, method = "BY"),
           prob_nonneg_adj = p.adjust( prob_nonneg, method = "BY"),
           prob_nonpos_adj = p.adjust( prob_nonpos, method = "BY")) %>%
    ungroup() %>%
    mutate(is_hit = p_value_adj <= 1E-2 & !(protgroup_id %in% exclude_protgroup_ids))
  fit_degr_contrast_stats[[df_name]] <- df
}

fit_rate_params.wide_stats <- reshape(fit_stats$rate_params %>%
                                      dplyr::select(protgroup_id, rate, param, condition, median = `50%`),
                                      direction="wide", idvar=c("protgroup_id", "rate", "param"),
                                      timevar="condition", v.names="median")

protgroup_info.df <- ms_data$protgroups %>%
  dplyr::inner_join(dplyr::select(fit_status.df, protgroup_id, is_converged = mcmc_converged)) %>%
  dplyr::mutate(long_label = paste0("Gene: ", gene_names, "<br>Protein: ", protein_names, "<br>ACs: ", majority_protein_acs, "<br>PG_ID#", protgroup_id),
                short_label = str_trunc(if_else(!is.na(gene_names), gene_names, majority_protein_acs), 20),
                is_relevant = !is_contaminant & !is_reverse & n_quanted > 0)

rate_cums.df <- dplyr::filter(fit_stats$rate_params, param == 'cum') %>%
  dplyr::mutate(median = `50%`) %>%
  dplyr::inner_join(dplyr::select(protgroup_info.df, protgroup_id, short_label, long_label, is_relevant)) %>%
  dplyr::mutate(median_log10 = log10(median),
                mean_log10 = log10(mean))
rate_cums_wide.df <- reshape(rate_cums.df %>% dplyr::select(condition, protgroup_id, rate, mean, mean_log10, median, median_log10, sd),
                             direction="wide", idvar=c("rate", "protgroup_id"),
                             timevar = "condition", v.names = c("mean", "median", "mean_log10", "median_log10", "sd")) %>%
  dplyr::inner_join(dplyr::select(protgroup_info.df, protgroup_id, short_label, long_label, is_relevant))

contrastXcondition.df <- dplyr::filter(as.data.frame(as.table(contrastXcondition)), Freq!=0.0) %>%
  dplyr::rename(weight = Freq) %>%
  dplyr::mutate(cond_role = if_else(weight > 0, "lhs", "rhs")) %>%
  dplyr::rename(contrast = contrasts, condition = conditions)

rate_cum_contrast_wide.df <- dplyr::group_by(contrastXcondition.df, contrast) %>% do({
  contr.df <- .
  contr_rate_cums.df <- dplyr::inner_join(rate_cums.df, contr.df)
  reshape(dplyr::select(contr_rate_cums.df, protgroup_id, is_relevant, short_label, cond_role, rate, mean, mean_log10, median, median_log10, sd),
          direction="wide", idvar=c("rate", "protgroup_id", "is_relevant", "short_label"),
          timevar = "cond_role", v.names = c("mean", "median", "mean_log10", "median_log10", "sd"))
}) %>% dplyr::ungroup() %>%
  dplyr::mutate(median_log10.lhs = log10(median.lhs),
                median_log10.rhs = log10(median.rhs))

# for renormalizing rate cums, skip otherwise
require(quantreg)

rate_cum_scales_wide.df <- dplyr::filter(protgroup_info.df, is_converged & n_quanted > 0.75 * max(n_quanted) & !is_contaminant & !is_reverse) %>%
  dplyr::select(protgroup_id) %>%
  dplyr::inner_join(rate_cums.df) %>%
  dplyr::group_by(rate) %>%
  do({
    df <- dplyr::mutate(., protgroup_id = factor(protgroup_id))
    glm_res <- rq(log10(median) ~ protgroup_id + condition, data=df,
                  contrast = list(condition="contr.treatment"))
    data.frame(intercept = glm_res$coefficients[["(Intercept)"]],
               shift.Mock = 0.0,
               shift.SC35M = glm_res$coefficients[["conditionSC35M"]],
               shift.SC35MdelNS1 = glm_res$coefficients[["conditionSC35MdelNS1"]],
               rho = glm_res$rho)
  }) %>% dplyr::ungroup()

rate_cum_scales.df <- reshape(rate_cum_scales_wide.df, direction="long",
                              idvar=c("rate"), timevar="condition", v.names="shift") %>%
  dplyr::mutate(mult = 10^shift)

rate_cum_contrast_stats.df <- dplyr::inner_join(rate_cum_scales.df, contrastXcondition.df) %>%
  dplyr::group_by(rate, contrast) %>%
  do({
    df <- dplyr::select(., contrast, cond_role, shift)
    dfw <- reshape(df, direction="wide", idvar="contrast", timevar="cond_role", v.names="shift")
    data.frame(intercept = 0.0,#glm_res$coefficients[["(Intercept)"]],
               slope_log10 = dfw$shift.rhs - dfw$shift.lhs,
               slope = 10^(dfw$shift.rhs - dfw$shift.lhs),
               stringsAsFactors=FALSE)
  })

message('Saving per-condition rate scales...')
save(project_id, job_version, rate_cum_scales.df,
     file = file.path(analysis_path, paste0(project_id, '_', job_version, '_param_scales.RData')))
# end of renormalizing
load(file.path(analysis_path, paste0(project_id, '_', job_version, '_param_scales_rescaled.RData')))

plots_path <- file.path(data_path, 'plots')

ggplot(dplyr::inner_join(fit_stats$recyc,
                         dplyr::select(ms_data$protgroups, protgroup_id, gene_names))) +
  geom_density(aes(x=`50%`, fill=mstag), alpha=0.5)

recyc_wide.df <- reshape(dplyr::select(fit_stats$recyc, protgroup_id, mstag, `25%`, `75%`, `50%`, mean),
                         direction="wide", v.names = c("mean", "25%", "50%", "75%"),
                         timevar = "mstag", idvar = "protgroup_id") %>%
  dplyr::inner_join(dplyr::select(ms_data$protgroups, majority_protein_acs, protgroup_id, gene_names, is_contaminant, is_reverse) %>%
                      dplyr::mutate(gene_label = str_trunc(if_else(!is.na(gene_names), gene_names, majority_protein_acs), 20),
                                    is_relevant = !is_contaminant & !is_reverse)) %>%
  dplyr::mutate(median_L.log10 = log10(`50%.L`),
                median_H.log10 = log10(`50%.H`))

recyc_kde2d <- kde2d4plot(recyc_wide.df, "median_L.log10", "median_H.log10", n=300)
recyc_wide.df <- recyc_kde2d$data

pdf(file=file.path(data_path, paste0("plots/", project_id, "_", job_version, "_recyc_rates.pdf")), width=13, height=12)
ggplot(recyc_wide.df, aes(x=`median_L.log10`, y=`median_H.log10`)) +
  geom_raster(data=dplyr::filter(recyc_kde2d$density_df, bin2d_density>0.01),
              aes(fill=bin2d_density^0.25)) +
  geom_contour(data=recyc_kde2d$density_df, stat="contour",
               aes(z=bin2d_density, color=..level..),
               breaks=c(0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 8.0, 10.0, 15.0, 20.0, 25.0, 30.0)) +
  geom_abline(data = data.frame(slope=1, intercept=0),
              aes(slope=slope, intercept=intercept),
              linetype=2, color="gray") +
  geom_point(data=dplyr::filter(recyc_wide.df, bin2d_density < 0.1), alpha=0.5) +
  geom_text(data=dplyr::filter(recyc_wide.df, bin2d_density < 0.1),
            aes(label=gene_label), color="firebrick", vjust=-1.1) +
  scale_color_gradient(low="gray75", high="black", guide=FALSE) +
  scale_fill_gradient("density", low="gray95", high="black") +
  scale_shape_manual(values=c("TRUE"=16, "FALSE"=1)) +
  theme_bw_ast(base_family="")
dev.off()

summary(glm(`50%.H` ~ `50%.L`, data=recyc_wide.df))

condition_palette <- c(SC35M="red", SC35MdelNS1="orange", Mock="darkgray")

pdf(file=file.path(data_path,  paste0("plots/", project_id, "_", job_version, "_recyc_rates_vs_synth.pdf")), width=14, height=12)
rates_vs_synth.df <- dplyr::filter(fit_stats$rate_params, var=="s_cum") %>%
  dplyr::rename(s_cum.median = `50%`) %>% dplyr::group_by(protgroup_id) %>%
  dplyr::filter(row_number(desc(s_cum.median)) == 1L) %>% # highest synthesis
  dplyr::inner_join(recyc_wide.df) %>%
  dplyr::inner_join(dplyr::filter(fit_stats$global, var=="L0_labu") %>%
                     dplyr::select(protgroup_id, `L0_labu.median` = `50%`))
ggplot(rates_vs_synth.df, aes(x=s_cum.median, y=pmin(0.1, `50%.L`))) +
  geom_point(aes(color=condition))+#color=pmax(-2, L0_labu.median))) +
  geom_text(data=dplyr::filter(rates_vs_synth.df, `50%.L` > 0.025 | s_cum.median > 1E+3),
            aes(label=gene_names), color="black", vjust=-1.1) +
  scale_x_log10() +
  #scale_color_continuous("L0.labu") +
  scale_color_manual(values=condition_palette) +
  theme_bw_ast(base_family="")
dev.off()

  ggplot(dplyr::filter(fit_stats$rate_params, var=="s_cum" & condition=="SC35MdelNS1") %>%
           dplyr::inner_join(recyc_wide.df) %>%
           dplyr::inner_join(dplyr::filter(fit_stats$global, var=="L0_labu") %>%
                               dplyr::select(protgroup_id, `50%.L0_labu` = `50%`)),
         aes(x=`50%.L0_labu`, y=`50%`)) +
    geom_point() +
    geom_text(data=dplyr::filter(fit_stats$rate_params, var=="s_cum" & condition=="SC35MdelNS1") %>%
                dplyr::inner_join(recyc_wide.df) %>%
                dplyr::inner_join(dplyr::filter(fit_stats$global, var=="L0_labu") %>%
                                    dplyr::select(protgroup_id, `50%.L0_labu` = `50%`)) %>%
                dplyr::filter(abs(log(`50%`) - `50%.L0_labu`) > 1.5),
              aes(label=gene_names), color="black", size=3, vjust=-1.1) +
    scale_y_log10()#
  facet_grid(. ~ condition)
  
rate_params_wide.df <- reshape(dplyr::select(fit_stats$rate_params, protgroup_id, condition, `25%`, `75%`, `50%`, mean, var),
                         direction="wide", v.names = c("mean", "25%", "50%", "75%"),
                         timevar = "var", idvar = c("protgroup_id", "condition")) %>%
  dplyr::inner_join(dplyr::select(ms_data$protgroups, protgroup_id, gene_names))

rate_params_wide_cond.df <- reshape(dplyr::select(fit_stats$rate_params, protgroup_id, condition, `25%`, `75%`, `50%`, mean, var),
                               direction="wide", v.names = c("mean", "25%", "50%", "75%"),
                               timevar = "condition", idvar = c("protgroup_id", "var")) %>%
  dplyr::inner_join(dplyr::select(ms_data$protgroups, protgroup_id, gene_names))

ggplot(rate_params_wide.df,
       aes(x=`50%.d0_cum`, y=`50%.d1_cum`)) +
  geom_point(alpha=0.1) +
  geom_text(data=dplyr::filter(rate_params_wide.df, abs(log10(`50%.d0_cum`)-log10(`50%.d1_cum`)) > 0.75),
            aes(label=gene_names), color="black", size=3, vjust=-1.1) +
  scale_x_log10() + scale_y_log10() +
  facet_grid(. ~ condition)

ggplot(dplyr::filter(rate_params_wide_cond.df, var %in% c('d0_cum', 'd1_cum')),
       aes(x=`50%.Mock`, y=`50%.SC35M`)) +
  geom_point() +
  geom_text(data=dplyr::filter(rate_params_wide_cond.df, var %in% c('d0_cum', 'd1_cum') &
                                 abs(log10(`50%.SC35M`)-log10(`50%.Mock`)) > 0.5),
            aes(label=gene_names), color="black", size=3, vjust=-1.1) +
  scale_x_log10() + scale_y_log10() +
  facet_grid(. ~ var)

sel_recyc_wide.df <- dplyr::filter(recyc_wide.df, `50%.H` > 0.02 | `50%.L` > 0.02)

ggplot(fit_cond_contrast_stats$rate_params %>% dplyr::filter(var %in% c('d0', 'd1') & param=='cum') %>%
       dplyr::inner_join(dplyr::select(ms_data$protgroups, protgroup_id, gene_names, majority_protein_acs)),
       aes(x=pmax(-2, pmin(2,`50%`)),
           y=pmin(prob_nonneg, prob_nonpos),
           label=str_trunc(if_else(!is.na(gene_names), gene_names, majority_protein_acs), 20))) +
  geom_point() +
  geom_text(data=fit_cond_contrast_stats$rate_params %>%
              dplyr::inner_join(dplyr::select(ms_data$protgroups, protgroup_id, gene_names, majority_protein_acs)) %>%
              dplyr::filter(var %in% c('d0', 'd1') & param=='cum' & pmin(prob_nonneg, prob_nonpos) <= 1E-5 & abs(`50%`) > 0.3),
            vjust=-1.1) +
  scale_y_continuous(trans = mlog10_trans()) +
  facet_wrap(~ contrast + var + param, scales="free") +
  theme_bw_ast(base_family = "")

protein_quant_stats.df <- pulse_data.df %>% group_by(protgroup_id) %>%
  dplyr::summarize(n_quant = sum(!is_imputed),
                   has_quant = any(!is_imputed))

ggplot(fit_rate_params.wide_stats %>% dplyr::inner_join(protein_quant_stats.df) %>%
       dplyr::filter(rate=='d0' & param == 'cum'),
       #aes(x=median.Mock/rate_scales.df$s_scale[1], y=median.SC35M/rate_scales.df$s_scale[2], color=n_quant)) +
       #aes(x=median.Mock/rate_scales.df$s_scale[1], y=median.SC35MdelNS1/rate_scales.df$s_scale[3], color=n_quant)) +
       aes(x=median.Mock, y=median.SC35M, color=n_quant)) +
  geom_point(alpha=0.1, size=1.0) +
  annotate('segment', x=0, y=0, xend=10, yend=10, linetype=2) +
  geom_text(data = fit_rate_params.wide_stats %>%
                   inner_join(proteins.df) %>%
                   dplyr::inner_join(protein_quant_stats.df) %>%
                   dplyr::filter(param=='cum' & rate=='d0'),
            aes(label=gene_names), size=3) +
  #scale_y_log10() + scale_x_log10() +
  xlim(0,1.0) + ylim(0,1.0) +
  facet_wrap(param ~ rate, scales = "free", ncol=3L)

fit_cond_contrast_stats$rate_params <- dplyr::mutate(
  fit_cond_contrast_stats$rate_params,
  p_value_fake = pmax(1E-250, rgamma(n(), shape=1E-2, scale=1E-125), p_value),
  p_value_fake2 = pmax(1E-100, rgamma(n(), shape=1E-2, scale=1E-50), p_value)
)
fit_degr_contrast_stats$rate_params <- dplyr::mutate(
  fit_degr_contrast_stats$rate_params,
  p_value_fake = pmax(1E-250, rgamma(n(), shape=1E-2, scale=1E-125), p_value),
  p_value_fake2 = pmax(1E-60, rgamma(n(), shape=1E-2, scale=1E-40), p_value)
)

pdf(file=file.path(data_path, "plots/d0_contrasts_volcanos.pdf"), width=12, height=18)
ggplot(fit_cond_contrast_stats$rate_params %>%
       dplyr::inner_join(contrast_quant_stats.df) %>%
       dplyr::filter(rate == 'd0' & param == 'cum'),
       aes(x=`50%`, y=p_value_fake, color=n_quant)) +
  geom_point(size=1, alpha=0.1) +
  geom_text(data = fit_cond_contrast_stats$rate_params %>%
              dplyr::inner_join(contrast_quant_stats.df) %>%
              dplyr::filter(param=='cum' & rate == 's' & p_value < 1E-10),
            aes(label=lead_gene_name), size=3,
            position = position_nudge(y=0.01)) +
  geom_vline(xintercept=0) +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_log10() +
  facet_grid(contrast ~ .)
dev.off()

pdf(file=file.path(data_path, "plots/d0_vs_d1_volcanos.pdf"), width=12, height=18)
ggplot(fit_degr_contrast_stats$rate_params %>% dplyr::filter(param=='cum'),
       aes(x=`50%`,
           #sign(`50%`)*log10(1.0+abs(`50%`)), 
           y=p_value_fake2, color=is_hit)) +
  geom_point() +
  geom_text(data = fit_degr_contrast_stats$rate_params %>%
                   dplyr::filter(is_hit & param=='cum'),
            aes(label=protgroup_id), vjust=-1.1, size=3) +
  geom_vline(xintercept=0) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="gray")) +
  scale_y_log10() +
  scale_x_continuous(limits = c(-10, 10)) +
  theme(legend.position="none")+
  facet_wrap(param ~ condition, scales = "free", ncol=1L)
dev.off()

write.table(pcp_protein_hits.df %>% dplyr::select(gene_names, description, protgroup_id, majority_protein_acs,
                                                  profile_subgroups, max_hits_rank, contrast, peak_p_value_adj.min,
                                                  is_single_peptide_subgroup ),
            file.path(analysis_path, paste0(job_dataset,'_',job_version,'_pcp_hits_20160606.txt')), sep='\t', row.names = F)

out_rdata_filepath <- file.path(analysis_path, paste0(project_id, '_', job_version, '_fit_all.RData'))
message('Saving full analysis results to ', out_rdata_filepath, '...')
analysis_info <- list(project_id = project_id, version = job_version)
save(analysis_info, contrastXcondition, contrastXcondition.df,
     fit_status.df, fit_cond_contrast_stats, fit_degr_contrast_stats, fit_stats,
     fit_rate_params.wide_stats,
     file = out_rdata_filepath)

require(rjson)
json_filepath <- file.path(analysis_path, paste0(project_id, '_', job_version, '_fit_stats.json'))
write(file = json_filepath, toJSON(fit_stats))

message( 'Dataset analysis has been saved' )
