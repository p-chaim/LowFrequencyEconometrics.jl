  push!(LOAD_PATH, "F:\\julia_files\\lfecon\\src")
  using lfecon
  using Distributions
  using PyPlot

  include("mw_lfst_pv.jl")
  include("mw_lrstat_vec.jl")
  include("mw_marg_lf.jl")
  include("mw_rw_cov.jl")
  include("mw_s_cov_ha.jl")
  include("mw_h_cov_ha_i0_i1_i2.jl") # ok
  include("mw_h_test_i0_i1_pv.jl") # ok
  #include("mw_stest_pv.jl")
  param_lfst_m = 10
  param_lfst_t = 20
  param_lfur_m = 14
  param_h = 6
  param_s = 5
  nsim = 100
  n_wsim = 100
  n_dsim = 100
  t1 = 500; t2 = 1000
  rnd_delta = randn(t1, t2)
  rmatrix_w_string = rnd_delta
  t1 = 1100; t2 = 1000
  rnd_w = randn(t1, t2)
  rmatrix_d_string =  rnd_w
  itrend = 1
  dat = Simul_ARMA(T = 123)
  x_t = psi_compute(length(dat), 12)'dat
  # Htest using I(0) Null
  ind_int = 0
  pv_h_i0, h_i0_stat=  mw_h_test_i0_i1_pv(x_t,param_h,itrend,ind_int,n_wsim,rmatrix_w_string,nsim)
  # Htest using UR Null @
  ind_int = 1
  pv_h_i1, h_i1_stat = mw_h_test_i0_i1_pv(x_t,param_h,itrend,ind_int,n_wsim,rmatrix_w_string,nsim)
# Stest
  q = size(x_t, 1)
  psi = psi_compute(100, q+1)#, itrend)
  psi = psi[:, 1:q]
  cov_h0_i0 = psi'psi
  cov_u = mw_rw_cov(100)
  cov_h0_i1 = psi'cov_u*psi

  pv_s_i0=mw_stest_pv(x_t,cov_h0_i0,param_s,n_dsim,rmatrix_d_string,nsim)
