#job.args <- c("phubel_pulse", "", "20170930", "0", "0", "447" )
if (!exists('job.args')) {
  job.args <- commandArgs(trailingOnly = TRUE)
}

project_id <- job.args[[1]]
message('Project ID=', project_id)

job_version <- job.args[[3]]
message('Job dataset version is ', job_version)

skip_existing <- job.args[[4]] > 0L
message("Skipping existing models: ", skip_existing)

job_id <- as.integer(job.args[[5]])
message('Job ID is ', job_id)

job_chunk <- as.integer(job.args[[6]])
sel_pg_id <- job_chunk - 1L
# job_chunk is output and processed later

source("~/R/config.R")
source(file.path(base_scripts_path, 'misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))
source(file.path(misc_scripts_path, 'reshape.R'))

skip_existing <- FALSE
is_stiff <- skip_existing

analysis_path <- scratch_path
out_rdata_path <- file.path(analysis_path, paste0("dynamics_fit_", job_version))
if (!dir.exists(out_rdata_path)) { dir.create(out_rdata_path) }
out_rdata_filepath <- file.path(out_rdata_path, paste0('pulse_dynamics_fit_', job_version, '_', sel_pg_id, '.RData'))
if (skip_existing) {
if (file.exists(out_rdata_filepath)) {
  stop("Stan model for ", sel_pg_id, " (", out_rdata_filepath, ") already exists, stopping")
} else {
  message("Stan model for ", sel_pg_id, " not found, starting inference...")
}
}

source(file.path(base_scripts_path, "msglm/R/stan_process_utils.R"))
source(file.path(base_scripts_path, "msglm/R/msglm_results.R"))
source(file.path(base_scripts_path, 'misc/ggplot_ext.R'))
source(file.path(misc_scripts_path, "julia_utils.R"))

mop.max_nprocesses <- 8
#mop.nprocesses <- 8
source(file.path(pipeline_scripts_path, 'init_cluster.R'))

require(dplyr)
require(rstan)
require(insilicoMop)
require(stringr)
require(ggplot2)

message("Load data...")
load(file.path(scratch_path, paste0(project_id, "_", job_version, "_lean.RData")))
load(file.path(scratch_path, "pulse_stan_models.RData"))

prepare_pulse_dynamics.stan_data <- function(pulse_data.df, t0=0.0, t_M_switch=12.0) {
  if (n_distinct(pulse_data.df$protgroup_id) > 1L) {
    stop("Data contains multiple protein groups: ", sort(unique(pulse_data.df$protgroup_id)))
  }
  timepoints = as.numeric(levels(ms_data$mschannels$timepoint))
  mschannels.df <- dplyr::select(ms_data$mschannels, mschannel, mstag, msrun, msrun_ix, condition, timepoint, replicate) %>%
    dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift=total_msrun_calib_shift))
  msruns.df <- dplyr::filter(mschannels.df, mstag=="Sum") %>%
    #dplyr::mutate(total_mschannel_shift = 0.0) %>%
    dplyr::arrange(msrun_ix)
  pulse_data_all.df <- expand.grid(protgroup_id = unique(pulse_data.df$protgroup_id),
                                   mschannel = unique(dplyr::filter(mschannels.df, mstag!="Sum")$mschannel)) %>%
    dplyr::left_join(pulse_data.df %>% dplyr::select(mschannel, protgroup_id, Intensity)) %>%
    dplyr::inner_join(dplyr::select(mschannels.df, mschannel, msrun, mstag, msrun_ix, condition, timepoint, replicate)) %>%
    dplyr::arrange(protgroup_id, msrun_ix, mstag) %>%
    #dplyr::mutate(Intensity = Intensity/exp(-msrun_shift)) %>%
    # remove M intensity: it should be zero, but because M switch function has to be continuous
    # it's not the case, so to avoid the side effects it just has to be ignored
    dplyr::filter(!(timepoint==timepoints[[1]] & mstag=="M")) %>%
    dplyr::mutate(observation_ix = row_number())
  if (all(is.na(pulse_data_all.df$Intensity))) {
    warning("Data contain no quantified intensities")
    median_intensity <- 1.0
  } else {
    median_intensity <- median(pulse_data.df$Intensity, na.rm = TRUE)
  }
  message('Median intensity=', median_intensity)
  if (t0 >= min(timepoints)) {
    warning("Integration starts (t0=", t0, ") show be strictly before the first data point (t_min=", min(timepoints), ")")
  }
  timepoints_sim <- as.numeric(1:42)
  if (n_distinct(pulse_data_all.df$protgroup_id) * n_distinct(pulse_data_all.df$condition) *
      n_distinct(pulse_data_all.df$timepoint) * n_distinct(pulse_data_all.df$replicate) *
      n_distinct(pulse_data_all.df$mstag) !=
      # should match obervations + skipped observations (M channel at M switch timepoint for all conditions)
      nrow(pulse_data_all.df) + n_distinct(pulse_data_all.df$condition)*n_distinct(pulse_data_all.df$replicate)) {
    stop("Data does not represent the expected 4D array")
  }
  list(
    Nrepl = n_distinct(pulse_data_all.df$replicate),
    Ncond = n_distinct(pulse_data_all.df$condition),
    Nt = length(timepoints),
    Nt_sim = sum(timepoints_sim >= t0),
    Nmsrun = n_distinct(pulse_data_all.df$msrun),
    t0 = t0,
    t_M_switch = t_M_switch,
    timepoints = timepoints,
    timepoints_sim = timepoints_sim[timepoints_sim >= t0],
    msrun_shift = msruns.df$total_msrun_shift,
    Nobservations = nrow(pulse_data_all.df),
    msrun2condition = as.integer(msruns.df$condition),
    msrun2timepoint = as.integer(msruns.df$timepoint), # timepoint should be a factor
    observation2msrun = as.array(pulse_data_all.df$msrun_ix),
    observation2label = as.array(as.integer(pulse_data_all.df$mstag)),
    Nquanted = sum(!is.na(pulse_data_all.df$Intensity)),
    quant2observation = as.array(dplyr::filter(pulse_data_all.df, !is.na(Intensity))$observation_ix),
    miss2observation = as.array(dplyr::filter(pulse_data_all.df, is.na((Intensity)))$observation_ix),
    q_data = as.array(dplyr::filter(pulse_data_all.df, !is.na(Intensity))$Intensity),
    global_labu_shift = global_labu_shift,
    g_params = growth_params.mtx,
    s_gscales = 1/rate_scales.df$s_scale,
    recyc = recycling_rates.df$rate
  )
}

sel_intensities.df <- dplyr::filter(ms_data$protgroup_intensities, protgroup_id == sel_pg_id)
                    #dplyr::filter(grepl("^TUBB$", gene_names, ignore.case = TRUE))
lead_protein_acs <- dplyr::semi_join(ms_data$protgroups, sel_intensities.df)$majority_protein_acs[[1]]
lead_protein_ac <- str_split_fixed(lead_protein_acs, fixed(';'), 2L)[[1]]
lead_gene_names <- dplyr::semi_join(ms_data$protgroups, sel_intensities.df)$gene_names[[1]]
lead_gene_name <- str_split_fixed(lead_gene_names, fixed(';'), 2L)[[1]]
message("Processing ", lead_gene_name, " (", lead_protein_ac, " PG=", sel_pg_id, ")")

#quantile(prot_pulse_data.df %>% dplyr::filter(mstag=='M' & timepoint > 12) %>% .$Intensity, 0.98) /
#  quantile(prot_pulse_data.df %>% dplyr::filter(mstag == 'L' & timepoint <= 12) %>% .$Intensity, 0.98)

pulse_dyn_stan_data <- prepare_pulse_dynamics.stan_data(sel_intensities.df, t0=0.0, t_M_switch = 12.0) %>%
  c(is_stiff = is_stiff) %>%
  c(def_norm_data)

# data.frame(q_l=log(pulse_dyn_stan_data$q_data),
#            obs_ix=pulse_dyn_stan_data$quant2observation,
#            msrun_ix=pulse_dyn_stan_data$observation2msrun[pulse_dyn_stan_data$quant2observation],
#            msrun_shift=pulse_dyn_stan_data$msrun_shift[pulse_dyn_stan_data$observation2msrun[pulse_dyn_stan_data$quant2observation]],
#            lbl_ix=pulse_dyn_stan_data$observation2label[pulse_dyn_stan_data$quant2observation],
#            cond_ix=pulse_dyn_stan_data$msrun2condition[pulse_dyn_stan_data$observation2msrun[pulse_dyn_stan_data$quant2observation]],
#            timept_ix=pulse_dyn_stan_data$msrun2timepoint[pulse_dyn_stan_data$observation2msrun[pulse_dyn_stan_data$quant2observation]]) %>%
#   View()

mcmc_start_time <- Sys.time()

message( 'Running STAN in HMC mode...' )
if ( !is.null( mop.cluster ) ) {
  stan_sampling_seed = sample.int(.Machine$integer.max, 1)
  clusterEvalQ( mop.cluster, library(rstan) )
  clusterEvalQ( mop.cluster, library(dplyr) )
  clusterExport( mop.cluster, varlist=c('stan_sampling_seed', 'pulse_dyn_stan_data', 'pulse_dynamics.stan_model') )
  pulse_dynamics.samples <- clusterApplyLB( mop.cluster, seq_len(mop.nprocesses), function( chain_ix ) {
    sampling( pulse_dynamics.stan_model,
              data = pulse_dyn_stan_data,
              #init = function() { pulse_dynamics.generate_init_params(pulse_dynamics.model_data) },
              iter=8000L, chain_id=chain_ix, seed=stan_sampling_seed, chains=8L, thin=4L )
  } )
  pulse_dynamics.samples <- sflist2stanfit( pulse_dynamics.samples )
} else {
  pulse_dynamics.samples <- sampling(pulse_dynamics.stan_model,
                                     data = pulse_dyn_stan_data,
                                     #init = function() { pulse_dynamics.generate_init_params(pulse_dynamics.model_data) },
                                     iter = 100L, chains = 4L, thin = 4L)
}

source( file.path( pipeline_scripts_path, 'stop_cluster.R' ) )

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
                               'd0_cum', 'd1_cum', 's_cum'),
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

min.iteration <- as.integer(1.0 * pulse_dynamics.samples@sim$warmup)

message('Computing parameters statistics...')
pulse_dynamics.stan_stats <- pulse_dynamics.samples %>%
  stan.extract_samples(pars = unlist(sapply(vars_info, function(vi) vi$names)),
                       min.iteration = min.iteration) %>%
  monitor( print = FALSE ) %>% as.data.frame
pulse_dynamics.stan_stats$var_name <- rownames(pulse_dynamics.stan_stats)

message('Extracting MCMC samples...')
pulse_dynamics.stan_samples <- stan.extract_samples(pulse_dynamics.samples,
                                                    pars = unlist(sapply(vars_info, function(vi) vi$names)),
                                                    min.iteration = min.iteration,
                                                    permuted = TRUE)

message('Preparing dimensions information...')
dim_info <- list( pos = data.frame(pos = factor(c('start', 'end'), levels=c('start', 'end'), ordered=TRUE), stringsAsFactors = FALSE ),
                  SILAC_LH = dplyr::select(ms_data$mschannels, mstag) %>% dplyr::filter(mstag %in% c('L','H')) %>% dplyr::distinct() %>% dplyr::arrange(mstag),
                  mstag = dplyr::select(ms_data$mschannels, mstag) %>% dplyr::distinct() %>% dplyr::arrange(mstag),
                  condition = dplyr::select(ms_data$mschannels, condition) %>% dplyr::distinct() %>% dplyr::arrange(condition),
                  time_sim = data.frame( timepoint = as.numeric(pulse_dyn_stan_data$timepoints_sim), stringsAsFactors = FALSE ) )

message('Composing model parameter reports...')
vars_results <- lapply(names(vars_info), vars_statistics, pulse_dynamics.stan_stats, pulse_dynamics.stan_samples, vars_info, dim_info)
names(vars_results) <- names(vars_info)
vars_results$rate_params$stats <- parse_rate_vars(vars_results$rate_params$stats)
vars_results$rate_points$stats <- parse_rate_vars(vars_results$rate_points$stats)

local({
  message( 'Calculating rate param condition-vs-condition contrasts' )
  samples.df <- reshape(dplyr::select(vars_results$rate_params$samples, -ends_with("bend"), -ends_with("steep")),
                        direction='long',
                        idvar=c('condition', 'chain', 'iteration', 'unpermuted_ix', 'sample_ix'),
                        timevar = 'param',
                        v.names = c('d0', 'd1', 's'), sep='_') #%>%
    #dplyr::left_join(rate_scales.df) %>%
    #dplyr::mutate(d0_scale = if_else(is.na(d0_scale), 1.0, d0_scale),
    #              d1_scale = if_else(is.na(d1_scale), 1.0, d1_scale),
    #              s_scale = if_else(is.na(s_scale), 1.0, s_scale),
    #              d0 = d0 / d0_scale,
    #              d1 = d1 / d1_scale,
    #              s = s / s_scale )
  contrast_stats.df <- vars_contrast_stats(samples.df, group_cols = 'param',
                                           var_names = c('d0', 'd1', 's'),
                                           condition_col = 'condition',
                                           contrastXcondition)
  #print(str(contrast_stats.df))
  vars_results$rate_params$condition_contrast_stats <<- contrast_stats.df %>% dplyr::mutate(rate = var)
  message( 'Calculating new-vs-old degradation rate contrasts' )
  d_samples.df <- reshape(samples.df %>% dplyr::rename(value.d0 = d0, value.d1=d1, value.s=s),
                          direction="long", idvar=c('condition', 'chain', 'iteration', 'unpermuted_ix', 'sample_ix', 'param'),
                          timevar="rate", v.names=c("value"), sep= ".") %>%
      dplyr::filter(rate %in% c('d0', 'd1'))

  d_contrast_stats.df <- vars_contrast_stats(d_samples.df, var_name = "value",
                                             group_cols = c('param', 'condition'),
                                             condition_col='rate',
                                             contrastXcondition = matrix(c(1, -1), dimnames=list(contrasts=c('d1_vs_d0'), rate=c('d1', 'd0')), ncol=2))
  #print(str(d_contrast_stats.df))
  vars_results$rate_params$degradation_contrast_stats <<- d_contrast_stats.df
})

message("Done calculating statistics")

# remove largest samples
vars_results$label_sim$samples <- NULL
vars_results$rate_sim$samples <- NULL

message( 'Preparing analysis report...' )
analysis_info = list(job_version = job_version, job_chunk = job_chunk,
                     protgroup_id = sel_pg_id,
                     gene_name = lead_gene_name,
                     protein_ac = lead_protein_ac)

# collect some convergence statistics
vars.df <- bind_rows(lapply(names(vars_results), function(frame_name) {
    data.frame(frame = frame_name,
               varname = unique(vars_results[[frame_name]]$stats$var),
               stringsAsFactors = FALSE)
}))
var_fit_status.df <- bind_rows(lapply(1:nrow(vars.df), function(var_ix) {
  dplyr::filter(vars_results[[vars.df$frame[[var_ix]] ]]$stats, var == vars.df$varname[[var_ix]]) %>%
    dplyr::group_by(var) %>%
    dplyr::summarise_at(.funs=funs(min=min(., na.rm=TRUE),
                               max=max(., na.rm=TRUE),
                               median=median(., na.rm=TRUE),
                               mean=mean(., na.rm=TRUE),
                               `25%`=quantile(., 0.25, na.rm=TRUE),
                               `75%`=quantile(., 0.75, na.rm=TRUE)),
                        .vars=vars(Rhat, n_eff)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(frame = vars.df$frame[[var_ix]])
}))

fit_status <- list(
  mcmc_converged = dplyr::filter(var_fit_status.df, var=="abu_sim" & frame=="label_sim")$Rhat_max <= 1.1,
  mcmc_start_time = mcmc_start_time,
  mcmc_end_time = Sys.time(),
  vars_status = var_fit_status.df
)
message("MCMC time: ", fit_status$mcmc_end_time - fit_status$mcmc_start_time)
print(fit_status$vars_status, n = -1)

message( 'Saving STAN results to ', out_rdata_filepath, '...' )
save(analysis_info, pulse_dyn_stan_data, vars_results, fit_status,
     file = out_rdata_filepath)
message( 'Done saving.' )

plots_path <- file.path(analysis_path, paste0("plots_", job_version))
plot_filename_prefix <- paste0(str_replace(lead_gene_name, "\\.", '_'), '_', lead_protein_ac, '_', sel_pg_id)
plot_title_prefix <- paste0(lead_gene_name, " (",  lead_protein_ac, ")" )

message("Generating plots...")

source(file.path(misc_scripts_path, 'ggplot_ext.R'))

condition_palette <- c(SC35M="red", SC35MdelNS1="orange", Mock="darkgray")

pdf(file.path(plots_path, paste0(plot_filename_prefix, "_data.pdf")), width=12, height=12)
shown_intensities.df <- dplyr::select(sel_intensities.df, mschannel, protgroup_id,
                                      Intensity, Intensity_msrun_norm, Intensity_total_calib_adj) %>%
  dplyr::left_join(ms_data$protgroups %>% dplyr::select(protgroup_id, gene_names)) %>%
  dplyr::left_join(ms_data$mschannels)
shown_intensities.df <- bind_rows(dplyr::mutate(shown_intensities.df, intensity=Intensity, norm_type="none"),
                                  dplyr::mutate(shown_intensities.df, intensity=Intensity_msrun_norm, norm_type="msrun"),
                                  dplyr::mutate(shown_intensities.df, intensity=Intensity_total_calib_adj, norm_type="total"))
ggplot(shown_intensities.df, aes(x=timepoint, y=intensity, color=condition, fill=condition)) +
  geom_point(aes(shape=factor(replicate))) +
  geom_smooth(aes(x=as.numeric(timepoint))) +
  scale_color_manual("condition", values=condition_palette) +
  scale_fill_manual("condition", values=condition_palette) +
  #scale_y_log10() +
  theme_bw_ast(base_family = "") +
  facet_grid(mstag + gene_names + protgroup_id ~ norm_type, scales = "free")
dev.off()

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

pdf(file.path(plots_path, paste0(plot_filename_prefix, "_dynamics.pdf")), width=10, height=8)
label_sim_stats_shown <- vars_results$label_sim$stats
label_sim_stats_max_shown <- label_sim_stats_shown %>%
  dplyr::filter(as.numeric(timepoint) >= 12.5) %>%
  dplyr::group_by(mstag) %>%
  dplyr::summarize(max_shown_sim = min(quantile(`75%`, 0.95) * 1.25, max(`97.5%`))) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(dplyr::filter(shown_intensities.df, norm_type=="total") %>%
                      dplyr::group_by(mstag) %>%
                      dplyr::summarise(max_shown_data = pmin(quantile(intensity, 0.9, na.rm=TRUE)*1.25, max(intensity, na.rm=TRUE))*exp(-global_labu_shift))) %>%
  dplyr::mutate(max_shown = pmax(max_shown_sim, pmin(2*max_shown_sim, max_shown_data, na.rm=TRUE), na.rm=TRUE))
label_sim_stats_shown <- mutate(label_sim_stats_shown %>% dplyr::left_join(label_sim_stats_max_shown),
                                `97.5%` = if_else(`97.5%` <= max_shown, `97.5%`, NA_real_ )) %>%
  mutate_at(funs(pmin(., max_shown)), .vars=vars(`2.5%`, `25%`, `50%`, `75%`)) %>%
  mutate_at(funs(Intensity=.*exp(global_labu_shift)), .vars=vars(`2.5%`, `25%`, `50%`, `75%`, `97.5%`))

ggplot(label_sim_stats_shown,
       aes(x=timepoint, fill=condition, color=condition)) +
  geom_ribbon(aes(ymin=`25%_Intensity`, ymax=`75%_Intensity`), alpha=0.5, stat = "identity" ) +
  geom_line(aes(y=`50%_Intensity`), size=1 ) +
  geom_line(aes(y=`2.5%_Intensity`), linetype=2 ) +
  geom_line(aes(y=`97.5%_Intensity`), linetype=2 ) +
  geom_vline(data=data.frame(timepoint=12), aes(xintercept=timepoint), linetype=2, color='gray25') +
  geom_point(data=dplyr::filter(shown_intensities.df, norm_type=="total") %>%
               dplyr::inner_join(label_sim_stats_max_shown) %>%
               dplyr::mutate(timepoint = as.numeric(as.character(timepoint))),
             aes(y=pmin(intensity, max_shown*exp(global_labu_shift))), alpha=0.25, size=0.50) +
  #annotate('vline', xintercept = 12, linetype=2, size=2, color='red') +
  scale_colour_manual('condition', values=condition_palette) +
  #scale_y_continuous('concentration',
  #    # starting amounts are poorly defined and jump a lot (esp. L label), so skip it
  #    limits=c(0, quantile(dplyr::filter(vars_results$label_sim$stats, timepoint >= 10)$`97.5%`, 0.9))) +
  scale_fill_manual('condition', values=condition_palette) +
  ggtitle(paste0(plot_title_prefix, " dynamics simulation")) +
  theme_bw_ast(base_family = "") +
  facet_grid(mstag ~ ., scales = "free")
dev.off()

rate_sim_stats_shown <- vars_results$rate_sim$stats
rate_sim_stats_max_shown <- rate_sim_stats_shown %>%
  dplyr::filter(as.numeric(timepoint) >= 12.5) %>%
  dplyr::group_by(var) %>%
  dplyr::summarize(max_shown = min(quantile(`75%`, 0.98) * 5, max(`97.5%`))) %>%
  dplyr::ungroup()
rate_sim_stats_shown <- mutate(rate_sim_stats_shown %>% dplyr::left_join(rate_sim_stats_max_shown),
                               `97.5%_limit` = if_else(`97.5%` <= max_shown, NA_real_, max_shown),
                               `97.5%` = if_else(`97.5%` <= max_shown, `97.5%`, NA_real_)) %>%
  mutate_at(funs(pmin(., max_shown)), .vars=vars(`2.5%`, `25%`, `50%`, `75%`))

pdf(file.path(plots_path, paste0(plot_filename_prefix, "_rate_dynamics.pdf")), width=10, height=8)
ggplot(rate_sim_stats_shown,
       aes(x=timepoint, fill=condition, color=condition)) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), alpha=0.5, stat = "identity") +
  geom_line(aes(y=`50%`), size=1) +
  geom_line(aes(y=`2.5%`), linetype=2) +
  geom_line(aes(y=`97.5%`), linetype=2) +
  geom_line(aes(y=`97.5%_limit`), linetype=3) +
  geom_vline(data=data.frame(timepoint=12), aes(xintercept=timepoint), linetype=2, color='gray25') +
  #annotate('vline', xintercept = 12, linetype=2, size=2, color='red') +
  scale_colour_manual('condition', values=condition_palette) +
  scale_fill_manual('condition', values=condition_palette) +
  ggtitle(paste0(plot_title_prefix, " rate dynamics")) +
  theme_bw_ast(base_family = "") +
  facet_grid(var ~ ., scales = "free")
dev.off()

d_rate_sim_stats_shown <- dplyr::filter(vars_results$rate_sim$stats, var %in% c('d0_sim', 'd1_sim'))
d_rate_sim_stats_max_shown <- d_rate_sim_stats_shown %>%
  dplyr::filter(as.numeric(timepoint) >= 12.5) %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(max_shown = min(quantile(`75%`, 0.98) * 5, max(`97.5%`))) %>%
  dplyr::ungroup()
d_rate_sim_stats_shown <- mutate(d_rate_sim_stats_shown %>% dplyr::left_join(d_rate_sim_stats_max_shown),
                                 `97.5%_limit` = if_else(`97.5%` <= max_shown, NA_real_, max_shown),
                                 `97.5%` = if_else(`97.5%` <= max_shown, `97.5%`, NA_real_)) %>%
  mutate_at(funs(pmin(., max_shown)), .vars=vars(`2.5%`, `25%`, `50%`, `75%`))

pdf(file.path(plots_path, paste0(plot_filename_prefix, "_d_rate_dynamics.pdf")), width=10, height=8)
ggplot(d_rate_sim_stats_shown %>% dplyr::filter(var %in% c('d0_sim', 'd1_sim')),
       aes(x=timepoint, color=var, fill=var)) +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), alpha=0.5, stat = "identity" ) +
  geom_line(aes(y=`50%`), size=1 ) +
  geom_line(aes(y=`2.5%`), linetype=2) +
  geom_line(aes(y=`97.5%`), linetype=2) +
  geom_line(aes(y=`97.5%_limit`), linetype=3) +
  geom_vline(data=data.frame(timepoint=12), aes(xintercept=timepoint), linetype=2, color='gray25') +
  #annotate('vline', xintercept = 12, linetype=2, size=2, color='red') +
  scale_colour_manual('rate', values=c('d0_sim'='darkgray','d1_sim'='skyblue3')) +
  scale_fill_manual('rate', values=c('d0_sim'='darkgray','d1_sim'='skyblue3')) +
  ggtitle(paste0(plot_title_prefix, " d0 vs d1 rate dynamics")) +
  theme_bw_ast(base_family = "") +
  facet_grid(condition ~ ., scales = "free")
dev.off()

message("Done plots")
