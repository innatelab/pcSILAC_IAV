project_id <- "phubel_pulse"
job_version <- "20171025"

source('~/R/config.R')
source(file.path(base_scripts_path, 'misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'reshape.R'))
source(file.path(base_scripts_path, 'misc/ggplot_ext.R'))

analysis_path <- scratch_path

require(ggplot2)
require(plotly)
require(dplyr)
require(stringr)

load(file.path(scratch_path, paste0(project_id, "_", job_version, "_lean.RData")))
load(file.path(analysis_path, paste0(project_id, '_', job_version, '_param_scales_rescaled.RData')))
load(file.path(analysis_path, paste0(project_id, '_', job_version, '_fit_all.RData')))

protgroup_info.df <- ms_data$protgroups %>%
  dplyr::inner_join(dplyr::select(fit_status.df, protgroup_id, is_converged = mcmc_converged)) %>%
  dplyr::mutate(long_label = paste0("Gene: ", gene_names, "<br>Protein: ", protein_names, "<br>ACs: ", majority_protein_acs, "<br>PG_ID#", protgroup_id),
                short_label = str_trunc(if_else(!is.na(gene_names), gene_names, majority_protein_acs), 20),
                is_relevant = !is_contaminant & !is_reverse)

# use cum_scaled after renormalizing
rate_cums.df <- dplyr::filter(fit_stats$rate_params, param == 'cum_scaled') %>%
  dplyr::mutate(median = `50%`) %>%
  dplyr::inner_join(dplyr::select(protgroup_info.df, protgroup_id, short_label, long_label, is_relevant)) %>%
  dplyr::mutate(median_log10 = log10(median),
                mean_log10 = log10(mean))
rate_cums_wide.df <- reshape(rate_cums.df %>% dplyr::select(condition, protgroup_id, rate, mean, mean_log10, median, median_log10, sd),
                             direction="wide", idvar=c("rate", "protgroup_id"),
                             timevar = "condition", v.names = c("mean", "median", "mean_log10", "median_log10", "sd")) %>%
  dplyr::inner_join(dplyr::select(protgroup_info.df, protgroup_id, short_label, long_label, is_relevant))

rate_cum_contrast_wide.df <- dplyr::group_by(contrastXcondition.df, contrast) %>% do({
  contr.df <- .
  contr_rate_cums.df <- dplyr::inner_join(rate_cums.df, contr.df)
  reshape(dplyr::select(contr_rate_cums.df, protgroup_id, is_relevant, short_label, cond_role, rate, mean, mean_log10, median, median_log10, sd),
          direction="wide", idvar=c("rate", "protgroup_id", "is_relevant", "short_label"),
          timevar = "cond_role", v.names = c("mean", "median", "mean_log10", "median_log10", "sd"))
}) %>% dplyr::ungroup() %>%
  dplyr::mutate(median_log10.lhs = log10(median.lhs),
                median_log10.rhs = log10(median.rhs))

# plot all rates contrasts
bin2d_cutoff <- 0.05
dplyr::mutate(rate_cum_contrast_wide.df, param="cum_scaled") %>%
dplyr::group_by(contrast, rate) %>% dplyr::do({
  sel_cums_wide.df <- .
  sel_rate <- sel_cums_wide.df$rate[[1]]
  sel_contr <- sel_cums_wide.df$contrast[[1]]
  message("Plotting ", sel_rate, " ", sel_contr, "...")
  cond_lhs <- dplyr::filter(contrastXcondition.df, contrast==sel_contr & cond_role=="lhs")$condition[[1]]
  cond_rhs <- dplyr::filter(contrastXcondition.df, contrast==sel_contr & cond_role=="rhs")$condition[[1]]
  sel_cums_wide_kde2d <- kde2d4plot(sel_cums_wide.df, "median_log10.lhs", "median_log10.rhs", n=300)
  sel_cums_wide.df <- dplyr::left_join(sel_cums_wide_kde2d$data,
                                       fit_cond_contrast_stats$rate_params) %>%
    dplyr::mutate(is_shown = p_value_adj < 0.005 & bin2d_density < bin2d_cutoff)
  #kde2d_breaks <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5)
  kde2d_breaks <- c(0.01, quantile(sel_cums_wide.df$bin2d_density[sel_cums_wide.df$bin2d_density>1E-5], c(seq(0.05, 0.95, by=0.05), 0.99)))
  kde2d_breaks <- sort(kde2d_breaks[kde2d_breaks >= 0.001 & kde2d_breaks < max(sel_cums_wide_kde2d$density_df$bin2d_density, na.rm=TRUE)])
  print(kde2d_breaks)
  pdf(file.path(data_path, "plots", paste0(project_id, "_", job_version, "_", sel_rate, "_", sel_contr, ".pdf")), width=12, height=10)
  print(ggplot(sel_cums_wide.df, aes(x=median_log10.rhs, y=median_log10.lhs)) +
    geom_raster(data=dplyr::filter(sel_cums_wide_kde2d$density_df, bin2d_density > 0.01),
                aes(fill=bin2d_density^0.25)) +
    geom_contour(data=sel_cums_wide_kde2d$density_df, stat="contour",
                 aes(z=bin2d_density, color=..level..),
                 breaks=kde2d_breaks) +
    geom_contour(data=sel_cums_wide_kde2d$density_df, stat="contour",
                 aes(z=bin2d_density),
                 breaks=bin2d_cutoff, color="red") +
    geom_point(data=dplyr::filter(sel_cums_wide.df, is_shown),
               aes(shape=is_relevant, size=-log10(if_else(is.na(p_value), 1E-200, pmax(p_value, 1E-200)))), alpha=0.5) +
    geom_text(data=dplyr::filter(sel_cums_wide.df, is_relevant & is_shown),
              aes(label=short_label), color="firebrick", vjust=-1.05, size=1.8) +
    geom_abline(data = data.frame(slope=1.0, intercept=0.0),
                aes(slope=slope, intercept=intercept), linetype=2, color="gray") +
    scale_color_gradient(low="gray75", high="black", guide=FALSE) +
    scale_fill_gradient("density", low="gray95", high="black") +
    scale_x_continuous(paste0("log10(median[",cond_rhs,"])")) +
    scale_y_continuous(paste0("log10(median[",cond_lhs,"])")) +
    scale_shape_manual(values=c("TRUE"=16, "FALSE"=1), guide=FALSE) +
    scale_size_continuous(range=c(2, 4), guide=FALSE) +
    theme_bw_ast(base_family="", base_size = 10))
  dev.off()
  data.frame()
})

# d1 vs d0 contrasts per condition
d_rate_cum_contrast_wide.df <- reshape(dplyr::select(rate_cums.df, protgroup_id, is_relevant, short_label, condition, rate, mean, mean_log10, median, median_log10, sd),
          direction="wide", idvar=c("protgroup_id", "condition", "is_relevant", "short_label"),
          timevar = "rate", v.names = c("mean", "median", "mean_log10", "median_log10", "sd")) %>%
  dplyr::mutate(median_log10.lhs = log10(median.d1),
                median_log10.rhs = log10(median.d0))

bin2d_cutoff <- 0.05
dplyr::mutate(d_rate_cum_contrast_wide.df, param="cum_scaled") %>%
  dplyr::group_by(condition) %>% dplyr::do({
    sel_cums_wide.df <- .
    sel_cond <- sel_cums_wide.df$condition[[1]]
    message("Plotting d0_vs_d1 ", sel_cond, "...")
    sel_cums_wide_kde2d <- kde2d4plot(sel_cums_wide.df, "median_log10.lhs", "median_log10.rhs", n=300)
    sel_cums_wide.df <- dplyr::left_join(sel_cums_wide_kde2d$data,
                                         fit_degr_contrast_stats$rate_params) %>%
      dplyr::mutate(is_shown = p_value_adj < 0.005 & bin2d_density < bin2d_cutoff)
    #kde2d_breaks <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5)
    kde2d_breaks <- c(0.01, quantile(sel_cums_wide.df$bin2d_density[sel_cums_wide.df$bin2d_density>1E-5], c(seq(0.05, 0.95, by=0.05), 0.99)))
    kde2d_breaks <- sort(kde2d_breaks[kde2d_breaks >= 0.001 & kde2d_breaks < max(sel_cums_wide_kde2d$density_df$bin2d_density, na.rm=TRUE)])
    print(kde2d_breaks)
    pdf(file.path(data_path, "plots", paste0(project_id, "_", job_version, "_d1_vs_d0_", sel_cond, ".pdf")), width=12, height=10)
    print(ggplot(sel_cums_wide.df, aes(x=median_log10.rhs, y=median_log10.lhs)) +
            geom_raster(data=dplyr::filter(sel_cums_wide_kde2d$density_df, bin2d_density > 0.01),
                        aes(fill=bin2d_density^0.25)) +
            geom_contour(data=sel_cums_wide_kde2d$density_df, stat="contour",
                         aes(z=bin2d_density, color=..level..),
                         breaks=kde2d_breaks) +
            geom_contour(data=sel_cums_wide_kde2d$density_df, stat="contour",
                         aes(z=bin2d_density),
                         breaks=bin2d_cutoff, color="red") +
            geom_point(data=dplyr::filter(sel_cums_wide.df, is_shown),
                       aes(shape=is_relevant, size=-log10(if_else(is.na(p_value), 1E-200, pmax(p_value, 1E-200)))), alpha=0.5) +
            geom_text(data=dplyr::filter(sel_cums_wide.df, is_relevant & is_shown),
                      aes(label=short_label), color="firebrick", vjust=-1.05, size=1.8) +
            geom_abline(data = data.frame(slope=1.0, intercept=0.0),
                        aes(slope=slope, intercept=intercept), linetype=2, color="gray") +
            scale_color_gradient(low="gray75", high="black", guide=FALSE) +
            scale_fill_gradient("density", low="gray95", high="black") +
            scale_x_continuous(paste0("log10(median[d0])")) +
            scale_y_continuous(paste0("log10(median[d1])")) +
            scale_shape_manual(values=c("TRUE"=16, "FALSE"=1), guide=FALSE) +
            scale_size_continuous(range=c(2, 4), guide=FALSE) +
            theme_bw_ast(base_family="", base_size = 10))
    dev.off()
    data.frame()
  })

rate_cum_contrasts.df <- dplyr::filter(fit_cond_contrast_stats$rate_params, param %in% c('cum', 'cum_scaled')) %>%
  dplyr::rename(contrast = contrasts) %>%
  dplyr::inner_join(dplyr::select(ms_data$protgroups, protgroup_id, is_contaminant, is_reverse)) %>%
  dplyr::mutate(orig_rate = rate,
                median = `50%`,
                p_value = 2*pmin(prob_nonneg, prob_nonpos),
                is_hit = p_value <= 0.001 & !is_contaminant & !is_reverse)
rate_cum_contrasts.df <- bind_rows(rate_cum_contrasts.df,
                                   dplyr::filter(rate_cum_contrasts.df, rate %in% c("d0", "d1")) %>% dplyr::group_by(protgroup_id, contrast, param) %>%
                                   dplyr::filter(row_number(-abs(median)) == 1L) %>% dplyr::ungroup() %>% dplyr::mutate(rate = "d"))
rate_cum_contrast_stats2.df <- dplyr::group_by(rate_cum_contrasts.df, rate, param) %>%
  dplyr::filter(is_hit) %>%
  dplyr::summarise(max_median = quantile(abs(median), 0.7)[[1]]) %>%
  dplyr::ungroup()
rate_cum_contrasts.df <- dplyr::inner_join(rate_cum_contrasts.df, rate_cum_contrast_stats2.df) %>%
  dplyr::mutate(median_clamped = pmin(pmax(-max_median, median), max_median))

rate_cum_contrasts_wide.df <- reshape(rate_cum_contrasts.df %>% dplyr::select(contrast, protgroup_id, p_value, is_hit, mean, rate, param, median, median_clamped, sd),
                                   direction="wide", idvar=c("contrast", "protgroup_id", "param"),
                                   timevar = "rate", v.names = c("p_value", "is_hit", "mean", "median", "median_clamped", "sd")) %>%
  #dplyr::filter(contrast != 'SC35M_vs_SC35MdelNS1') %>%
  dplyr::mutate(d_shown = if_else(abs(median.d0) == abs(median.d), "d0", "d1")) %>%
  dplyr::inner_join(dplyr::select(protgroup_info.df, protgroup_id, gene_names, protein_names, majority_protein_acs, is_contaminant, is_reverse, is_viral, is_relevant)) %>%
  dplyr::mutate(hit_type = str_replace_all(paste(if_else(is_hit.s, if_else(median.s > 0, 'Synth↑', 'Synth↓'), ''),
                                           if_else(is_hit.d, if_else(median.d > 0, 'Degr.↑', 'Degr.↓'), '')),
                                           "^\\s+|\\s+$", ""),
                long_label = paste0("Gene: ", gene_names, "<br>Protein: ", protein_names, "<br>ACs: ", majority_protein_acs, "<br>PG_ID#", protgroup_id),
                short_label = str_trunc(if_else(!is.na(gene_names), gene_names, majority_protein_acs), 20)
  ) %>%
  dplyr::mutate(hit_type = if_else(hit_type != "", hit_type, "neutral"))

hit_type_levels <- c("Synth↑", "Synth↑ Degr.↓", "Degr.↓", "Synth↓ Degr.↓",
                     "Synth↓", "Synth↓ Degr.↑", "Degr.↑", "Synth↑ Degr.↑")

pie_rate_cum_contrasts_wide.df <- dplyr::mutate(rate_cum_contrasts_wide.df,
                                             pos_hit_type = factor(if_else(hit_type == 'neutral', sample(hit_type_levels, n(), replace=TRUE), hit_type),
                                                                   levels=hit_type_levels, ordered=TRUE),
                                             hit_type = factor(if_else(hit_type=="", "neutral", hit_type),
                                                               levels=c(hit_type_levels, "neutral"), ordered=TRUE)) %>%
  dplyr::filter(hit_type != "neutral" | runif(n()) < 0.001)

hit_type_colors = c("Synth↑" = "gold", "Synth↑ Degr.↑" = "orange",
                    "Degr.↑" = "orangered", "Synth↓ Degr.↑" = "plum",
                    "Synth↓" = "green", "Synth↓ Degr.↓" = "lightskyblue",
                    "Degr.↓" = "blue", "Synth↑ Degr.↓" = "greenyellow",
                    "neutral" = "gray")
cairo_pdf(file.path(data_path, "plots", paste0(project_id, "_", job_version, "_S_vs_D_windrose.pdf")), width=9, height=6)
ggplot(pie_rate_cum_contrasts_wide.df %>% dplyr::filter(!is_viral & hit_type != "neutral")) +
  geom_bar(aes(x=pos_hit_type, fill=hit_type), width=1.0) +
  geom_text(data = pie_rate_cum_contrasts_wide.df %>% dplyr::filter(!is_viral & hit_type != "neutral") %>%
              dplyr::group_by(pos_hit_type, contrast, hit_type, param) %>% dplyr::summarise(n_protgroups = n()) %>% dplyr::ungroup(),
            aes(x=pos_hit_type, y=pmax(2, n_protgroups * 0.1), label=n_protgroups), color="white") +
  scale_x_discrete(drop=FALSE) +
  #scale_size_manual(values = c("TRUE"=1.0,"FALSE"=5)) +
  scale_fill_manual(values=hit_type_colors) +
  #scale_color_manual(values=hit_type_colors) +
  coord_polar(start = pi/2-pi/8, direction=1) +
  facet_grid(param ~ contrast, scales = "free") +
  theme_bw_ast(base_family = "Helvetica", base_size = 10) +
  scale_y_log10() +
  theme(legend.position="none")
dev.off()

#ggplotly(
pdf(file.path(data_path, "plots", paste0(project_id, "_", job_version, "_S_vs_D1.pdf")), width = 18, height=10)
ggplot(rate_cum_contrasts_wide.df %>% dplyr::filter(!is_hit.s & !is_hit.d),
       aes(x=median_clamped.s, y=median_clamped.d)) +
  geom_point(color="gray", alpha=0.5, size=1, shape=1) +
  geom_point(data=rate_cum_contrasts_wide.df %>% dplyr::filter(is_hit.s | is_hit.d),
             aes(color = hit_type, label=short_label, shape=is_viral)) +
  #geom_text(data=rate_cum_contrasts_wide.df %>%
  #              #dplyr::filter(str_detect(gene_names, "\\b(SQSTM1|CLUH|GPKOW|LGALS3BP)\\b"))
  #              dplyr::filter(is_hit.s & abs(median.s) > 0.1 | is_hit.d & abs(median.d) > 0.1)
  #          , aes(label=short_label, color=hit_type), vjust=-1.1, size=2) +
  scale_x_continuous() + scale_y_continuous() +
  scale_color_manual(values=hit_type_colors)+
  scale_shape_manual(values=c("TRUE"=17, "FALSE"=16))+
  facet_grid(. ~ contrast) +
  theme_bw_ast(base_family = "", base_size = 10)
dev.off()
#)

rate_cum_contrasts_wide_kde2d <- kde2d4plot(rate_cum_contrasts_wide.df, "median_clamped.s", "median_clamped.d", n=30)
rate_cum_contrasts_wide.df <- rate_cum_contrasts_wide_kde2d$data %>%
      dplyr::mutate(is_shown = pmin(p_value.s, p_value.d1) < 0.005 & bin2d_density < 0.001)
#kde2d_breaks <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5)
kde2d_breaks <- c(0.01, quantile(rate_cum_contrasts_wide.df$bin2d_density[!is.na(rate_cum_contrasts_wide.df$bin2d_density) & rate_cum_contrasts_wide.df$bin2d_density>1E-15], c(seq(0.05, 0.95, by=0.05), 0.99)))
kde2d_breaks <- sort(kde2d_breaks[kde2d_breaks >= 0.001 & kde2d_breaks < max(rate_cum_contrasts_wide_kde2d$density_df$bin2d_density, na.rm=TRUE)])
print(kde2d_breaks)
pdf(file.path(data_path, "plots", paste0(project_id, "_", job_version, "_d1_vs_d0_", sel_cond, ".pdf")), width=12, height=10)
ggplot(rate_cum_contrasts_wide.df,
             aes(x=median_clamped.s, y=median_clamped.d)) +
        geom_raster(data=dplyr::filter(rate_cum_contrasts_wide_kde2d$density_df, bin2d_density > 0.01),
                    aes(fill=bin2d_density^0.25)) +
        #geom_contour(data=rate_cum_contrasts_wide_kde2d$density_df, stat="contour",
        #             aes(z=bin2d_density, color=..level..),
        #             breaks=kde2d_breaks) +
        #geom_contour(data=rate_cum_contrasts_wide_kde2d$density_df, stat="contour",
        #             aes(z=bin2d_density),
        #             breaks=0.001, color="red") +
        #geom_point(data=dplyr::filter(rate_cum_contrasts_wide.df, is_shown),
        #           aes(shape=is_relevant, size=-log10(if_else(is.na(pmin(p_value.s, p_value.d)), 1E-200, pmax(pmin(p_value.s, p_value.d), 1E-200)))), alpha=0.5) +
        #geom_text(data=dplyr::filter(rate_cum_contrasts_wide.df, is_relevant & is_shown),
        #          aes(label=short_label), color="firebrick", vjust=-1.05, size=1.8) +
        scale_color_gradient(low="gray75", high="black", guide=FALSE) +
        scale_fill_gradient("density", low="gray95", high="black") +
        scale_x_continuous(paste0("log10(median[d0])")) +
        scale_y_continuous(paste0("log10(median[d1])")) +
        scale_shape_manual(values=c("TRUE"=16, "FALSE"=1), guide=FALSE) +
        scale_size_continuous(range=c(2, 4), guide=FALSE) +
        theme_bw_ast(base_family="", base_size = 10)
dev.off()


