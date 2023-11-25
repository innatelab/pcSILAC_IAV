# Build Stan models for PCP data
# 
# Author: Alexey Stukalov
###############################################################################

source('~/R/config.R')
source(file.path(base_scripts_path, 'misc/setup_base_paths.R'))

require(rstan)
pulse_dynamics.stan_model <- stan_model(file=file.path(base_scripts_path, 'adhoc/phubel_pulse/pulse_silac_dynamics.stan'),
                                        model_name = 'pulse_silac_dynamic', save_dso=TRUE, auto_write=TRUE)

rdata_filepath <- file.path(base_scratch_path, "phubel_pulse", paste0('pulse_stan_models.RData'))
message('Saving Stan models to ', rdata_filepath, '...')
save(pulse_dynamics.stan_model, file = rdata_filepath)
