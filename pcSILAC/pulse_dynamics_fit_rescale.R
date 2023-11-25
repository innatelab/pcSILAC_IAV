#job.args <- c("phubel_pulse", "", "20170930", "10", "3364" )
if (!exists('job.args')) {
  job.args <- commandArgs(trailingOnly = TRUE)
}

project_id <- job.args[[1]]
message('Project ID=', project_id)

job_version <- job.args[[3]]
message('Job dataset version is ', job_version)

job_id <- as.integer(job.args[[4]])
message('Job ID is ', job_id)

job_chunk <- as.integer(job.args[[5]])
sel_pg_id <- job_chunk - 1L
# job_chunk is output and processed later

source("~/R/config.R")
source(file.path(base_scripts_path, 'misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))
source(file.path(misc_scripts_path, 'reshape.R'))

analysis_path <- scratch_path
in_path <- file.path(analysis_path, paste0("dynamics_fit_", job_version))
in_rdata_filepath <- file.path(in_path, paste0('pulse_dynamics_fit_', job_version, '_', sel_pg_id, '.RData'))
if (!file.exists(in_rdata_filepath)) {
  stop("Stan model for ", sel_pg_id, " (", in_rdata_filepath, ") not found, stopping")
} else {
  message("Stan model for ", sel_pg_id, " found, starting fixing...")
}
out_path <- file.path(analysis_path, paste0("dynamics_fit_rescaled_", job_version))
if (!dir.exists(out_path)) dir.create(out_path)
out_rdata_filepath <- file.path(out_path, paste0('pulse_dynamics_fit_', job_version, '_', sel_pg_id, '.RData'))

source(file.path(base_scripts_path, "msglm/R/stan_process_utils.R"))
source(file.path(base_scripts_path, "msglm/R/msglm_results.R"))
source(file.path(base_scripts_path, 'misc/ggplot_ext.R'))
source(file.path(misc_scripts_path, "julia_utils.R"))

require(dplyr)
require(rstan)
require(insilicoMop)
require(stringr)
require(ggplot2)
#require(rjson)

message("Load data...")
load(file.path(scratch_path, paste0(project_id, "_", job_version, "_lean.RData")))

message("Load rate scales...")
load(file.path(scratch_path, paste0(project_id, '_', job_version, '_param_scales.RData')))

sel_intensities.df <- dplyr::filter(ms_data$protgroup_intensities, protgroup_id == sel_pg_id)
#dplyr::filter(grepl("^TUBB$", gene_names, ignore.case = TRUE))
lead_protein_acs <- dplyr::semi_join(ms_data$protgroups, sel_intensities.df)$majority_protein_acs[[1]]
lead_protein_ac <- str_split_fixed(lead_protein_acs, fixed(';'), 2L)[[1]]
lead_gene_names <- dplyr::semi_join(ms_data$protgroups, sel_intensities.df)$gene_names[[1]]
lead_gene_name <- str_split_fixed(lead_gene_names, fixed(';'), 2L)[[1]]
message("Processing ", lead_gene_name, " (", lead_protein_ac, " PG=", sel_pg_id, ")")

message("Load model fit...")
load(in_rdata_filepath)

timepoints_sim <- sort(unique(vars_results$rate_sim$stats$timepoint))

sigmoid <- function(t, p1, p2, p3, p4) {
  return(p3/(1+exp(-(p1 + t*p2))) + p4)
}

sigmoid_params <- function(args, start_val, end_val) {
  fa <- 1/(1+exp(args[[1]]))
  fb <- 1/(1+exp(-args[[2]]))

  c(-args[[1]], args[[2]]+args[[1]],
    (end_val-start_val)/(fb-fa),
    start_val-(end_val-start_val)/(fb-fa)*fa)
}

sigmoid_integral_params <- function(args, start_val, end_val) {
  a0 <- -args[[1]]
  a1 <- args[[2]]
  exp0p1 <- exp(a0)+1.0
  exp1p1 <- exp(-a1)+1.0
  expdm1 <- 1-exp(a0-a1)
  k <- (end_val-start_val)/(a1-a0)*exp0p1*exp1p1/expdm1
  
  c(sigmoid_params(args, start_val, end_val),
    (start_val*exp0p1 - end_val*exp(a0)*exp1p1)/expdm1,
    k, -k*log1p(exp(a0)))
}

sigmoid_integral <- function(t, p) {
  p[[5]]*t + p[[6]]*log1p(exp(p[[2]]*t+p[[1]])) + p[[7]]
}

rate_cum <- function(start_arg, end_arg, start_val, end_val) {
  a0 <- -start_arg
  a1 <- end_arg
  exp0p1 <- exp(a0)+1.0
  exp1p1 <- exp(-a1)+1.0
  # analytic form of \int_{0}^{1} smoother()dt
  ((start_val*exp0p1 - end_val*exp(a0)*exp1p1) +
    (log1p(exp(a1))-log1p(exp(a0)))*(end_val-start_val)/(a1-a0)*exp0p1*exp1p1)/(1.0-exp(a0-a1))
}

rate_points_samples.df <- vars_results$rate_points$samples
rate_points_samples_wide.df <- reshape(rate_points_samples.df, direction="wide",
                                       idvar = c("sample_ix", "chain", "iteration", "unpermuted_ix", "index_condition"),
                                       v.names = c("d0_args", "d1_args", "s_args", "d0_vals", "d1_vals", "s_vals"),
                                       drop=c("pos"),
                                       timevar="index_pos")

rate_cum_scales.df <- dplyr::mutate(rate_cum_scales.df, inv_mult = 10^(-shift))

cum_samples_wide.df <- dplyr::inner_join(vars_results$rate_params$samples, rate_points_samples_wide.df) %>%
  dplyr::inner_join(reshape(dplyr::select(rate_cum_scales.df, condition, rate, inv_mult),
                            direction="wide", idvar="condition", timevar = "rate", v.names = "inv_mult")) %>%
  dplyr::mutate(d0_cum = rate_cum(d0_args.1, d0_args.2, d0_vals.1, d0_vals.2)*max(timepoints_sim),
                d1_cum = rate_cum(d1_args.1, d1_args.2, d1_vals.1, d1_vals.2)*max(timepoints_sim),
                s_cum = rate_cum(s_args.1, s_args.2, s_vals.1, s_vals.2)*max(timepoints_sim),
                d0_cum_scaled = d0_cum*inv_mult.d0,
                d1_cum_scaled = d1_cum*inv_mult.d1,
                s_cum_scaled = s_cum*inv_mult.s)

# add scaled cum
vars_results$rate_params$samples <- dplyr::select(vars_results$rate_params$samples, -ends_with("cum_scaled")) %>%
  dplyr::inner_join(dplyr::select(cum_samples_wide.df, unpermuted_ix, index_condition, ends_with("cum_scaled")))

vars_info <- list(
  global = list(names = c('L0_labu',
                          'd_arg_sigma', 's_arg_sigma',
                          'd0_arg_delta_tau', 'd1_arg_delta_tau', 's_arg_delta_tau',
                          'd_val_sigma', 's_val_sigma',
                          'd0_val_delta_tau', 'd0_end_val_delta_tau',
                          'd1_val_delta_tau',
                          's_val_delta_tau', 's_end_val_delta_tau',
                          'd0_end_val_delta_lambda', 's_end_val_delta_lambda'),
                dims = c()),
  rate_points = list(names = c('d0_args', 'd1_args', 's_args',
                               'd0_vals', 'd1_vals', 's_vals'),
                     dims = c('pos', 'condition')),
  rate_params = list(names = c('d0_bend', 'd1_bend', 's_bend',
                               'd0_steep', 'd1_steep', 's_steep',
                               'd0_cum', 'd1_cum', 's_cum',
                               'd0_cum_scaled', 'd1_cum_scaled', 's_cum_scaled'),
                     dims = c('condition')),
  #recyc = list(names = c('recyc'), dims = c('SILAC_LH')),
  label_sim = list(names = c('abu_sim'),
                   dims = c('condition', 'time_sim', 'mstag') ),
  rate_sim = list(names = c('d0_sim', 'd1_sim', 's_sim'),
                  dims = c('time_sim', 'condition'))
)

parse_rate_vars <- function(rate_df) {
  dplyr::mutate(rate_df,
                rate = str_match(var, "^(d0|d1|s)_")[,2],
                param = str_match(var,"^(d0|d1|s)_(.+)$")[,3],
                param = recode(param, args="arg", vals="val"))
}

message('Assembling rescaled samples...')
cum_samples.df <- reshape(dplyr::rename(cum_samples_wide.df,
                                        val.d0_cum_scaled=d0_cum_scaled,
                                        val.d1_cum_scaled=d1_cum_scaled,
                                        val.s_cum_scaled=s_cum_scaled),
                            direction="long", idvar=c("index_condition", "unpermuted_ix", "chain", "iteration"),
                            timevar="var", v.names = "val") %>%
    dplyr::mutate(var_name=paste0(var, "[", index_condition, "]")) %>%
    dplyr::select(var, var_name, index_condition, unpermuted_ix, chain, iteration, val) %>%
  dplyr::arrange(var, var_name, unpermuted_ix)

message('Computing parameters statistics...')
cum_samples.arr <- array(cum_samples.df$val,
                         dim=c(n_distinct(cum_samples.df$iteration), n_distinct(cum_samples.df$chain), n_distinct(cum_samples.df$var_name)),
                         dimnames = list(iterations=NULL, chains=NULL, parameters=unique(cum_samples.df$var_name)))
cum_scaled.stan_stats <- as.data.frame(monitor(cum_samples.arr, print=FALSE))
cum_scaled.stan_stats$var_name <- rownames(cum_scaled.stan_stats)

message('Extracting MCMC samples...')

cum_scaled.stan_samples <- lapply(unique(cum_samples.df$var), function(cur_var) {
  var_samples.df <- dplyr::filter(cum_samples.df, var==cur_var) %>%
      dplyr::arrange(index_condition, unpermuted_ix)
  array(var_samples.df$val,
        dim=c(n_distinct(var_samples.df$index_condition),
              n_distinct(var_samples.df$unpermuted_ix)))
})
names(cum_scaled.stan_samples) <- unique(cum_samples.df$var)

message('Preparing dimensions information...')
dim_info <- list( pos = data.frame(pos = factor(c('start', 'end'), levels=c('start', 'end'), ordered=TRUE), stringsAsFactors = FALSE ),
                  SILAC_LH = dplyr::select(ms_data$mschannels, mstag) %>% dplyr::filter(mstag %in% c('L','H')) %>% dplyr::distinct() %>% dplyr::arrange(mstag),
                  mstag = dplyr::select(ms_data$mschannels, mstag) %>% dplyr::distinct() %>% dplyr::arrange(mstag),
                  condition = dplyr::select(ms_data$mschannels, condition) %>% dplyr::distinct() %>% dplyr::arrange(condition),
                  time_sim = data.frame( timepoint = as.numeric(pulse_dyn_stan_data$timepoints_sim), stringsAsFactors = FALSE ) )

message('Composing model parameter reports...')
cum_scaled_vars_results <- lapply(names(vars_info), vars_statistics, cum_scaled.stan_stats, cum_scaled.stan_samples, vars_info, dim_info)
names(cum_scaled_vars_results) <- names(vars_info)
cum_scaled_vars_results$rate_params$stats <- parse_rate_vars(cum_scaled_vars_results$rate_params$stats)

vars_results$rate_params$stats <- dplyr::anti_join(vars_results$rate_params$stats,
                                                   dplyr::select(cum_scaled_vars_results$rate_params$stats, var_name)) %>%
  dplyr::bind_rows(., cum_scaled_vars_results$rate_params$stats)
# FIXME rate sim not supported
#vars_results$rate_sim$stats <- dplyr::bind_rows(dplyr::anti_join(vars_results$rate_sim$stats,
#                                                                 dplyr::select(cum_scaled_vars_results$rate_sim$stats, var_name)),
#                                                cum_scaled_vars_results$rate_sim$stats)

local({
  message( 'Calculating rate param condition-vs-condition contrasts' )
  samples.df <- reshape(dplyr::select(vars_results$rate_params$samples, -ends_with("bend"), -ends_with("steep")),
                        direction='long',
                        idvar=c('condition', 'chain', 'iteration', 'unpermuted_ix', 'sample_ix'),
                        timevar = 'param',
                        v.names = c('d0', 'd1', 's'), sep='_')
  contrast_stats.df <- vars_contrast_stats(samples.df, group_cols = 'param',
                                           var_names = c('d0', 'd1', 's'),
                                           condition_col = 'condition',
                                           contrastXcondition)
  #print(str(contrast_stats.df))
  cum_scaled_vars_results$rate_params$condition_contrast_stats <<- dplyr::ungroup(contrast_stats.df %>% dplyr::mutate(rate = var))
  message( 'Calculating new-vs-old degradation rate contrasts' )
  d_samples.df <- reshape(samples.df %>% dplyr::rename(value.d0 = d0, value.d1=d1, value.s=s),
                          direction="long", idvar=c('condition', 'chain', 'iteration', 'unpermuted_ix', 'sample_ix', 'param'),
                          timevar="rate", v.names=c("value"), sep= ".") %>%
    dplyr::filter(rate %in% c('d0', 'd1'))
  
  d_contrast_stats.df <- vars_contrast_stats(d_samples.df, var_name = "value",
                                             group_cols = c('param', 'condition'),
                                             condition_col='rate',
                                             contrastXcondition = matrix(c(1, -1), dimnames=list(contrast=c('d1_vs_d0'), rate=c('d1', 'd0')), ncol=2))
  #print(str(d_contrast_stats.df))
  cum_scaled_vars_results$rate_params$degradation_contrast_stats <<- dplyr::ungroup(d_contrast_stats.df)
})
vars_results$rate_params$condition_contrast_stats <- dplyr::ungroup(vars_results$rate_params$condition_contrast_stats)
vars_results$rate_params$condition_contrast_stats <- dplyr::bind_rows(
  dplyr::anti_join(vars_results$rate_params$condition_contrast_stats,
                   dplyr::select(cum_scaled_vars_results$rate_params$condition_contrast_stats, param, index_contrast, var, rate)),
  cum_scaled_vars_results$rate_params$condition_contrast_stats)

vars_results$rate_params$degradation_contrast_stats <- dplyr::ungroup(vars_results$rate_params$degradation_contrast_stats)
vars_results$rate_params$degradation_contrast_stats <- dplyr::bind_rows(
  dplyr::anti_join(vars_results$rate_params$degradation_contrast_stats,
                   dplyr::select(cum_scaled_vars_results$rate_params$degradation_contrast_stats, param, index_contrast, var)),
  cum_scaled_vars_results$rate_params$degradation_contrast_stats)

message("Done calculating statistics")

# remove largest samples
vars_results$label_sim$samples <- NULL
vars_results$rate_sim$samples <- NULL

message( 'Saving STAN results to ', out_rdata_filepath, '...' )
save(analysis_info, pulse_dyn_stan_data, vars_results, fit_status,
     file = out_rdata_filepath)
message('Done saving.')

plots_path <- file.path(analysis_path, paste0("plots_", job_version))
plot_filename_prefix <- paste0(str_replace(lead_gene_name, "\\.", '_'), '_', lead_protein_ac, '_', sel_pg_id)
plot_title_prefix <- paste0(lead_gene_name, " (",  lead_protein_ac, ")" )

message("Generating plots...")

source(file.path(misc_scripts_path, 'ggplot_ext.R'))

condition_palette <- c(SC35M="red", SC35MdelNS1="orange", Mock="darkgray")

pdf(file.path(plots_path, paste0(plot_filename_prefix, "_params.pdf")), width=10, height=8)
ggplot(vars_results$rate_params$stats %>% dplyr::mutate(var_label = paste0(param, "(", rate, ")")) %>%
       dplyr::group_by(rate, param) %>% dplyr::mutate(`97.5%_limit` = max(`50%` + (`75%`-`25%`)*5)) %>% dplyr::ungroup()) +
  geom_boxplot(aes(x=condition, middle=`50%`, lower=`25%`, upper=`75%`, ymin=`2.5%`, ymax=pmin(`97.5%`, `97.5%_limit`),
                   color=condition, fill=condition), alpha=0.5, stat = "identity" ) +
  scale_colour_manual('condition', values=condition_palette) +
  scale_fill_manual('condition', values=condition_palette) +
  facet_wrap(~ var_label, scales = "free", ncol=3) +
  theme_bw_ast(base_family = "") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  ggtitle(paste0(plot_title_prefix, " dynamics simulation"))
dev.off()

message("Done plots")
