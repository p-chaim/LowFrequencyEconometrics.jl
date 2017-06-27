function mw_stest_pv(x_t, cov_h0, param, n_dsim, rnd_delta, nsim)
  q = size(x_t, 1)
  cov_ha = mw_s_cov_ha(cov_h0, param, rnd_delta)
  # Construct Test Statistic
  lr = mw_lrstat_vec(x_t, cov_h0, cov_ha)
  xsim = randn(q, size(cov_ha, 3))
  a = chol(cov_h0)
  xsim = a'xsim
  lr_sim = mw_lrstat_vec(xsim, cov_h0, cov_ha)
  return lr, lr_sim
end
