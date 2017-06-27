function mw_lrstat_vec(x, cov, cov_a)
  nsim = size(x, 2)
  q = size(x, 1)
  # Compute likelihood under ho @
  lf_h0 = mw_marg_lf(x, cov)
  lf_ha = NaN*zeros(size(x,2), size(cov_a, 3))
  for i ∈ 1:size(cov_a, 3)
    lf_ha[:, i] = mw_marg_lf(x, cov_a[:,:,i])
  end
  lf_ha_avg = mean(lf_ha, 2)
  lrstat = lf_ha_avg./lf_h0
  return lrstat
end
