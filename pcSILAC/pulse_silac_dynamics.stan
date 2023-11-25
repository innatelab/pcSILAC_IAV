functions {

  real sigmoid(real t, real[] p) {
    real r;
    r = inv_logit(p[1]+t*p[2])*p[3] + p[4];
    return r;
  }

  real[] sigmoid_params(vector args, real start_val, real end_val) {
    real p[4];
    real fa;
    real fb;
    p[1] = -args[1];
    p[2] = args[2]+args[1];
    fa = inv_logit(-args[1]);
    fb = inv_logit(args[2]);
    p[3] = (end_val-start_val)/(fb-fa);
    p[4] = start_val-p[3]*fa;
    return p;
  }

  real[] sigmoid_integral_params(vector args, real start_val, real end_val)
  {
    real a0;
    real a1;
    real exp0p1;
    real exp1p1;
    real expdm1;
    real k;
    real p[7];

    a0 = -args[1];
    a1 = args[2];
    exp0p1= exp(a0)+1.0;
    exp1p1 = exp(-a1)+1.0;
    expdm1 = 1-exp(a0-a1);
    k = (end_val-start_val)/(a1-a0)*exp0p1*exp1p1/expdm1;

    p[1:4] = sigmoid_params(args, start_val, end_val);
    p[5] = (start_val*exp0p1 - end_val*exp(a0)*exp1p1)/expdm1;
    p[6] = k;
    p[7] = -k*log1p_exp(a0);
    return p;
  }

  real sigmoid_integral(real t, real[] p)
  {
      return p[5]*t + p[6]*log1p_exp(p[2]*t+p[1]) + p[7];
  }

  real softmax0(real x, real s) {
    return log_sum_exp(s*x, 0)/s;
  }

  real[] pulse_ode_rhs(real t, real[] Labels, real[] theta, real[] x_r, int[] x_i) {
    real d0_rate;
    real d1_rate;
    real s_rate;
    real g_rate;
    real g;
    real recycL;
    real recycH;
    real s_M_scale;
    real dLabels_dt[4];
    real t_max;
    real t_M_switch;
    int Nparams;
    int NgrowthParams;

    Nparams = x_i[1];
    NgrowthParams = x_i[2];
    t_max = x_r[1];
    t_M_switch = x_r[2];
    g_rate = sigmoid(t/t_max, theta[(3*Nparams+2+1):(3*Nparams+2+NgrowthParams)])/t_max;
    g = exp(sigmoid_integral(t/t_max, theta[(3*Nparams+2+1):(3*Nparams+2+NgrowthParams)]));
    d0_rate = sigmoid(t/t_max, theta[1:Nparams]);
    d1_rate = sigmoid(t/t_max, theta[(Nparams+1):(2*Nparams)]);
    s_rate = sigmoid(t/t_max, theta[(2*Nparams+1):(3*Nparams)]);
    recycL = theta[3*Nparams+1];
    recycH = theta[3*Nparams+2];
    s_M_scale = inv_logit(50.0*(t-t_M_switch));
    dLabels_dt[1] = -(d0_rate+g_rate)*Labels[1];
    dLabels_dt[2] = -(d1_rate+g_rate)*Labels[2] + recycL*s_rate;
    dLabels_dt[3] = -(d1_rate+g_rate)*Labels[3] + (s_M_scale*recycH+(1-s_M_scale))*s_rate;
    dLabels_dt[4] = -(d1_rate+g_rate)*Labels[4] + s_M_scale*s_rate;

    return dLabels_dt;
  }

  matrix pulse_dynamics(int is_stiff, real t0, real L0, real Lmin,
                         vector d0_args, real d0_start_val, real d0_end_val,
                         vector d1_args, real d1_start_val, real d1_end_val,
                         vector s_args, real s_start_val, real s_end_val,
                         vector recyc, real[] timepoints,
                         real t_M_switch,
                         real[] g_int_params, real[] x_r, int[] x_i)
  {
    int  Nt;
    int  Nparams;
    int  NgrowthParams;
    real theta[3*x_i[1]+2+x_i[2]];      // ODE parameters
    real xLabels0[4];            // initial state
    real xLabels[size(timepoints), 4];
    matrix[3, size(timepoints)] Labels;

    Nparams = x_i[1];
    NgrowthParams = x_i[2];
    Nt = size(timepoints);
    xLabels0[1] = L0;
    xLabels0[2] = 0;
    xLabels0[3] = 0;
    xLabels0[4] = 0;

    theta[1:Nparams] = sigmoid_params(d0_args, d0_start_val, d0_end_val);
    theta[(Nparams+1):(2*Nparams)] = sigmoid_params(d1_args, d1_start_val, d1_end_val);
    theta[(2*Nparams+1):(3*Nparams)] = sigmoid_params(s_args, s_start_val, s_end_val);
    theta[3*Nparams+1] = recyc[1];
    theta[3*Nparams+2] = recyc[2];
    theta[(3*Nparams+2+1):(3*Nparams+2+NgrowthParams)] = g_int_params;

    if (is_stiff == 0) {
        xLabels = integrate_ode_rk45(pulse_ode_rhs, xLabels0, t0, timepoints, theta, x_r, x_i, 1E-2, 1E-3, 2000);
    } else {
        xLabels = integrate_ode_bdf(pulse_ode_rhs, xLabels0, t0, timepoints, theta, x_r, x_i, 1E-2, 1E-3, 2000);
    }
    // xLabels[1] (L old) + xLabels[2] (L new)
    Labels[1,:] = Lmin + to_row_vector(xLabels[:,1]) + to_row_vector(xLabels[:,2]);
    Labels[2,:] = Lmin + to_row_vector(xLabels[:,3]);
    Labels[3,:] = Lmin + to_row_vector(xLabels[:,4]);
    
    return Labels;
  }

  real rate_cum(vector args, real start_val, real end_val)
  {
    real a0;
    real a1;
    real exp0p1;
    real exp1p1;

    a0 = -args[1];
    a1 = args[2];
    exp0p1 = exp(a0)+1.0;
    exp1p1 = exp(-a1)+1.0;
    // analytic form of \int_{0}^{1} smoother()dt
    return ((start_val*exp0p1 - end_val*exp(a0)*exp1p1) +
        (log1p_exp(a1)-log1p_exp(a0))*(end_val-start_val)/(a1-a0)*exp0p1*exp1p1)/(1.0-exp(a0-a1));
  }

  real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
      return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
  }
}

data {
  int<lower=0,upper=1> is_stiff;
  int<lower=1> Ncond;
  int<lower=1> Nt;
  int<lower=0> Nt_sim;
  int<lower=1> Nmsrun;
  real<lower=0> t0;         // when to start integrating
  real<lower=0> t_M_switch; // when medium switched to M
  real<lower=0> timepoints[Nt]; // timepoints
  real<lower=0> timepoints_sim[Nt_sim]; // timepoints for simulation
  matrix[4,Ncond] g_params;          // parameters for cell growth rate sigmoid
  vector<lower=0>[2] recyc;   // recycling coefficients for L and H

  vector[Nmsrun] msrun_shift;
  int<lower=0> Nobservations;
  int<lower=1,upper=Ncond> msrun2condition[Nmsrun];
  int<lower=1,upper=Nt> msrun2timepoint[Nmsrun];
  int<lower=1,upper=Nmsrun> observation2msrun[Nobservations];
  int<lower=1,upper=3> observation2label[Nobservations];

  int<lower=0,upper=Nobservations> Nquanted;
  int<lower=1,upper=Nobservations> quant2observation[Nquanted];
  int<lower=1,upper=Nobservations> miss2observation[Nobservations-Nquanted];
  real<lower=0> q_data[Nquanted]; // quantitations ordered by condition and msrun

  // global model parameters
  real global_labu_shift;

  // instrument calibrated parameters 
  real<lower=0> zDetectionFactor;
  real zDetectionIntercept;
  real<lower=0, upper=1> detectionMax;

  real<lower=0> sigmaScaleHi;
  real<lower=0> sigmaScaleLo;
  real sigmaOffset;
  real sigmaBend;
  real sigmaSmooth;

  real zShift;
  real zScale;
}

transformed data {
  real mzShift; // zShift for the missing observation intensity (zShift shifted by global_labu_shift)
  vector[size(q_data)] zScore; // log(q_data) transformed to z-score
  vector[size(q_data)] q_lsd; // log(sd(q_data))-global_labu_shift
  vector[size(q_data)] q_data_scaled; // q_data/exp(global_labu_shift)
  vector<lower=0>[size(q_data)] q_tau_scaled; // exp(msrun_shift+global_labu_shift)/sd(q_data)=exp(msrun_shift-q_lsd)
  vector<lower=0>[size(q_data)] q_scaled; // q_data/sd(q_data)
  vector[size(q_data)] q_shift; // msrun_shifts[]
  vector[Nobservations-size(q_data)] m_shift; // msrun_shifts[]
  real g_int_params[Ncond,7];          // parameters for cell growth rate sigmoid integrap

  int<lower=1,upper=Ncond*Nt*3> observation2prediction[Nobservations];
  int<lower=1,upper=Ncond*Nt*3> quant2prediction[Nquanted];
  int<lower=1,upper=Ncond*Nt*3> miss2prediction[Nobservations-Nquanted];
  int<lower=1,upper=Nt> observation2timepoint[Nobservations];

  real<lower=0> Lmin;           // minimal label intensity (to avoid bernoulli_logit() NaNs)
  real L0_labu_prior;
  real x_r[2];                  // real data for ODE system
  int x_i[2];                   // integer data for ODE system
  int<lower=2> NrateVars;       // the number of variables to model d0, d1 or s for a given condition
  int<lower=2> NrateParams;     // the number of coefficients for smoother(t,...) to define rate curve
  int<lower=2> NgrowthParams;   // the number of parameters for growth rate and growth integral

  real t_sim_max;
  real t_max;
  real qSum_t1;   // sum of all intensities at t1
  int Nquant_t1;  // N of quantitations at t1

  for (i in 1:Nobservations) {
    int msrun;
    msrun = observation2msrun[i];
    observation2prediction[i] = (msrun2condition[msrun]-1)*Nt*3 + (msrun2timepoint[msrun]-1)*3 + observation2label[i];
  }

  q_shift = msrun_shift[observation2msrun[quant2observation]];
  observation2timepoint = msrun2timepoint[observation2msrun];
  {
    vector[size(q_data)] q_ldata;
    int lastCond;

    q_ldata = log(to_vector(q_data));
    zScore = (q_ldata - zShift) * zScale;
    mzShift = zShift - global_labu_shift;

    // process the intensity data to optimize likelihood calculation
    for (i in 1:size(q_data)) {
      q_lsd[i] = intensity_log_std(zScore[i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
    }
    q_data_scaled = to_vector(q_data) ./ exp(global_labu_shift);
    q_tau_scaled = exp(q_shift + global_labu_shift - q_lsd);
    for (i in 1:size(q_data)) {
      if (q_tau_scaled[i] > 20.0) {
        q_tau_scaled[i] = 20.0;
      }
    }
    q_scaled = exp(q_ldata - q_lsd);

    lastCond = 0; // check that msruns are sorted by condition
    for (i in 1:Nmsrun) {
      int cond;
      cond = msrun2condition[i];
      if (cond < lastCond) {
        print("condition2msrun[",i,"]=",cond, ", while the last condition was #", lastCond);
      }
      lastCond = cond;
    }

    // calculate intensity sums for T1
    qSum_t1 = 0;
    Nquant_t1 = 0;
    for (i in 1:size(q_data)) {
      int obs;
      int msrun;
      obs = quant2observation[i];
      msrun = observation2msrun[obs];
      if (msrun2timepoint[msrun] == 1) {
        qSum_t1 = qSum_t1 + exp(q_ldata[i] - msrun_shift[msrun] - global_labu_shift);
        Nquant_t1 = Nquant_t1 + 1;
      }
    }
  }
  m_shift = msrun_shift[observation2msrun[miss2observation]];
  quant2prediction = observation2prediction[quant2observation];
  miss2prediction = observation2prediction[miss2observation];
  print("quant2prediction=", quant2prediction);
  print("miss2prediction=", miss2prediction);

  t_max = max(timepoints);     // T max
  t_sim_max = max(timepoints_sim);
  x_r[1] = t_max;
  x_r[2] = t_M_switch;
  for (i in 1:Ncond) {
    g_int_params[i, :] = sigmoid_integral_params(sub_col(g_params, 1, i, 2), g_params[3, i], g_params[4, i]);
    print("g[",i,"](0)=", exp(sigmoid_integral(0.0, g_int_params[i, :])));
    print("g[",i,"](Tmax)=", exp(sigmoid_integral(1.0, g_int_params[i, :])));
  }

  NrateVars = 4;
  NrateParams = 4;
  NgrowthParams = 7;
  x_i[1] = NrateParams;
  x_i[2] = NgrowthParams;
  // make sure L0 contains something very small nonzero to degeneration of bernoulli_logit() at t0
  L0_labu_prior = log(2*(qSum_t1 > 0 ? qSum_t1 : 1E-3)/max(1, Nquant_t1)); // back-extrapolate L0 from t1 to t0, 2=Nlabels at t1
  Lmin = exp(-10.0-zShift+global_labu_shift);
  print("q_lsd=", q_lsd);
  print("q_tau_scaled=", q_tau_scaled);
  print("q_scaled=", q_scaled);
  print("L0=", exp(L0_labu_prior));
  print("L[min]=", Lmin);
  print("data timepoints=", timepoints);
  print("sim timepoints=", timepoints_sim);
  print("observation2timepoint=", observation2timepoint);
}

parameters {
  matrix<lower=0.1,upper=10>[2, Ncond] d0_args; // condition-specific degradation rate start/end segment (after infection)
  matrix<lower=0.1,upper=10>[2, Ncond] d1_args; // condition-specific degradation rate start/end segment (before infection)
  matrix<lower=0.1,upper=10>[2, Ncond] s_args;  // condition-specific synthesis rate start/end segment

  real<lower=0.0> d_start_val;
  real<lower=0.0> s_start_val;
  vector<lower=0.0>[Ncond] d0_end_vals;
  vector<lower=0.0>[Ncond] d1_end_vals;
  vector<lower=0.0>[Ncond] s_end_vals;

  //vector<lower=0,upper=0.25>[2] recyc;  // recycling coefficients
  real L0_labu;              // starting amount of L label
  real<lower=0.25> d_arg_sigma;
  real<lower=0.25> s_arg_sigma;

  real<lower=0> d0_arg_delta_tau;
  matrix<lower=0>[2,Ncond-1] d0_arg_delta_lambda;
  real<lower=0> d1_arg_delta_tau;
  vector<lower=0>[2*Ncond] d1_arg_delta_lambda;
  real<lower=0> s_arg_delta_tau;
  matrix<lower=0>[2, Ncond-1] s_arg_delta_lambda;

  real<lower=0> d_val_sigma;
  real<lower=0> s_val_sigma;
  real<lower=0.01> d0_val_delta_tau;
  // 0.01>0 to avoid degenerated solutions with tau and lambda close to zero
  vector<lower=0.01>[Ncond-1] d0_val_delta_lambda;
  real<lower=0.01> d0_end_val_delta_tau;
  real<lower=0.01> d0_end_val_delta_lambda;
  real<lower=0.01> d1_val_delta_tau;
  vector<lower=0.01>[Ncond] d1_val_delta_lambda;
  real<lower=0.01> s_val_delta_tau;
  vector<lower=0.01>[Ncond-1] s_val_delta_lambda;
  real<lower=0.01> s_end_val_delta_tau;
  real<lower=0.01> s_end_val_delta_lambda;
}

transformed parameters {
}

model {
  // base condition rates
  col(d0_args, 1) ~ normal(0, d_arg_sigma);
  col(d1_args, 1) ~ normal(0, d_arg_sigma);
  col(s_args, 1) ~ normal(0, s_arg_sigma);

  // d0 and s for all conditions like the first condition
  d_start_val ~ normal(0, d_val_sigma);
  d0_end_vals ~ normal(0, d_val_sigma);
  d0_end_vals[1] - d_start_val ~ normal(0, d0_end_val_delta_tau*d0_end_val_delta_lambda);
  d0_end_vals[2:Ncond] - d0_end_vals[1] ~ normal(0, d0_val_delta_tau*d0_val_delta_lambda);

  s_start_val ~ normal(0, s_val_sigma);
  s_end_vals ~ normal(0, s_val_sigma);
  s_end_vals[1] - s_start_val ~ normal(0, s_end_val_delta_tau*s_end_val_delta_lambda);
  s_end_vals[2:Ncond] - s_end_vals[1] ~ normal(0, s_val_delta_tau*s_val_delta_lambda);

  for (i in 1:2) {
    d0_args[i,2:Ncond] - d0_args[i,1] ~ normal(0.0, d0_arg_delta_tau*d0_arg_delta_lambda[i,:]);
    s_args[i,2:Ncond] - s_args[i,1] ~ normal(0.0, s_arg_delta_tau*s_arg_delta_lambda[i,:]);
  }

  // d1 rate like d0 rate
  d1_end_vals ~ normal(0, d_val_sigma);
  d1_end_vals - d0_end_vals ~ normal(0.0, d1_val_delta_tau*d1_val_delta_lambda);
  to_vector(d1_args) - to_vector(d0_args) ~ normal(0.0, d1_arg_delta_tau*d1_arg_delta_lambda);

  L0_labu ~ cauchy(L0_labu_prior, 10.0);
  //recyc ~ normal(0, 0.01);

  d_arg_sigma ~ cauchy(0, 1.0);
  s_arg_sigma ~ cauchy(0, 1.0);
  d0_arg_delta_tau ~ cauchy(0, 1.0);
  to_vector(d0_arg_delta_lambda) ~ cauchy(0, 1.0);
  d1_arg_delta_tau ~ cauchy(0, 1.0);
  d1_arg_delta_lambda ~ cauchy(0, 1.0);
  s_arg_delta_tau ~ cauchy(0, 1.0);
  to_vector(s_arg_delta_lambda) ~ cauchy(0, 1.0);

  d_val_sigma ~ cauchy(0, 0.01);
  s_val_sigma ~ cauchy(0, 0.1);
  d0_val_delta_tau ~ cauchy(0, 1.0);
  d0_val_delta_lambda ~ cauchy(0, 1.0);
  d0_end_val_delta_tau ~ cauchy(0, 1.0);
  d0_end_val_delta_lambda ~ cauchy(0, 1.0);
  d1_val_delta_tau ~ cauchy(0, 0.1);
  d1_val_delta_lambda ~ cauchy(0, 1.0);
  s_val_delta_tau ~ cauchy(0, 1.0);
  s_val_delta_lambda ~ cauchy(0, 1.0);
  s_end_val_delta_tau ~ cauchy(0, 1.0);
  s_end_val_delta_lambda ~ cauchy(0, 1.0);

  //print("d0_args=",d0_args);
  {
    vector[Ncond*Nt*3] abu_hat_vec;
    vector[size(quant2observation)] q_pred_scaled;
    vector[size(quant2observation)] q_labu;
    vector[size(miss2observation)] m_labu;

    // model the dynamics in each condition and fill corresponding q_labu and m_labu blocks
    for (k in 1:Ncond) {
      matrix[3,Nt] abu_hat;// predicted normalized abundance for a given condition
      abu_hat = pulse_dynamics(is_stiff, t0, exp(L0_labu), Lmin,
                               col(d0_args, k), d_start_val, d0_end_vals[k],
                               col(d1_args, k), d_start_val, d1_end_vals[k],
                               col(s_args, k), s_start_val, s_end_vals[k],
                               recyc, timepoints, t_M_switch, g_int_params[k, :], x_r, x_i);
      abu_hat_vec[1+Nt*3*(k-1):Nt*3*k] = to_vector(abu_hat);
      //print("recyc[L]=", recyc[1], " recyc[H]=", recyc[2]);
      //print("labu_hat[",k,"]=", labu_hat);
    }
    q_pred_scaled = abu_hat_vec[quant2prediction].*q_tau_scaled;
    q_labu = q_shift + log(abu_hat_vec[quant2prediction]);
    m_labu = m_shift + log(abu_hat_vec[miss2prediction]);

    //print("L0=", exp(L0_labu));
    //print("abu_hat_vec=", abu_hat_vec);
    //print("q_pred_scaled=", q_pred_scaled);
    q_scaled ~ double_exponential(q_pred_scaled, 1.0);
    // model non-missing data
    1 ~ bernoulli_logit(q_labu * (zScale * zDetectionFactor) - mzShift * zScale * zDetectionFactor + zDetectionIntercept);
    // model missing data
    0 ~ bernoulli_logit(m_labu * (zScale * zDetectionFactor) - mzShift * zScale * zDetectionFactor + zDetectionIntercept);
  }
}

generated quantities {
  matrix<lower=0>[2,Ncond] d0_vals;
  matrix<lower=0>[2,Ncond] d1_vals;
  matrix<lower=0>[2,Ncond] s_vals;
  real abu_sim[Ncond, Nt_sim, 4];
  real d0_sim[Nt_sim, Ncond];
  real d1_sim[Nt_sim, Ncond];
  real s_sim[Nt_sim, Ncond];
  real d0_cum[Ncond];
  real d1_cum[Ncond];
  real s_cum[Ncond];
  real d0_steep[Ncond];
  real d1_steep[Ncond];
  real s_steep[Ncond];
  real d0_bend[Ncond];
  real d1_bend[Ncond];
  real s_bend[Ncond];

  for (k in 1:Ncond) {
    real d0_params[NrateParams];
    real d1_params[NrateParams];
    real s_params[NrateParams];

    d0_vals[1,k] = d_start_val;
    d0_vals[2,k] = d0_end_vals[k];
    d1_vals[1,k] = d_start_val;
    d1_vals[2,k] = d1_end_vals[k];
    s_vals[1,k] = s_start_val;
    s_vals[2,k] = s_end_vals[k];
    d0_params = sigmoid_params(col(d0_args, k), d0_vals[1,k], d0_vals[2,k]);
    d1_params = sigmoid_params(col(d1_args, k), d1_vals[1,k], d1_vals[2,k]);
    s_params = sigmoid_params(col(s_args, k), s_vals[1,k], s_vals[2,k]);

    abu_sim[k,:,1:3] = to_array_2d(pulse_dynamics(is_stiff, t0, exp(L0_labu), Lmin,
                                          col(d0_args, k), d0_vals[1,k], d0_vals[2,k],
                                          col(d1_args, k), d1_vals[1,k], d1_vals[2,k],
                                          col(s_args, k), s_vals[1,k], s_vals[2,k],
                                          recyc, timepoints_sim, t_M_switch, g_int_params[k, :], x_r, x_i)');
    for (i in 1:Nt_sim) {
      real t;
      t = (timepoints_sim[i]-0)/(timepoints_sim[Nt_sim]-0);
      d0_sim[i, k] = sigmoid(t, d0_params);
      d1_sim[i, k] = sigmoid(t, d1_params);
      s_sim[i, k] = sigmoid(t, s_params);
      abu_sim[k,i,4] = sum(abu_sim[k,i,1:3]);
    }
    d0_cum[k] = rate_cum(col(d0_args, k), d0_vals[1,k], d0_vals[2,k])*t_max;
    d1_cum[k] = rate_cum(col(d1_args, k), d1_vals[1,k], d1_vals[2,k])*t_max;
    s_cum[k] = rate_cum(col(s_args, k), s_vals[1,k], s_vals[2,k])*t_max;
    d0_steep[k] = d0_args[1,k] + d0_args[2,k];
    d1_steep[k] = d1_args[1,k] + d1_args[2,k];
    s_steep[k] = s_args[1,k] + s_args[2,k];
    d0_bend[k] = d0_args[1,k]/d0_steep[k] - 0.5;
    d1_bend[k] = d1_args[1,k]/d1_steep[k] - 0.5;
    s_bend[k] = s_args[1,k]/s_steep[k] - 0.5;
  }
}
