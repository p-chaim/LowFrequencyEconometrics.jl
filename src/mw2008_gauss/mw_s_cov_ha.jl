function mw_s_cov_ha(cov_h0, param, rnd_delta; big=true)#, rmatrix_string)
  # s_cov_ha - Generates nsim HA covariance matrices for S-test
  # param = local parameter
  # nsim = number of simulations
  # rmatrix_string = name of matrix with standard random variables (on disk)
  nsim = size(rnd_delta, 2)
  q = size(cov_h0, 1)
  param = 5.0
  qq =Int( q*(q+1)/2 )
  delta = rnd_delta
  delta = cumsum(delta[1:q, 1:nsim])/sqrt(q)
  delta = delta*(param/sqrt(q))
  gam = exp.(delta)
  cov_ha = NaN*zeros(qq, nsim)
  big_cov_ha = NaN*zeros(q, q, nsim)
  for i ∈ 1:nsim
    tmp = (cov_h0.*gam[:, i]).*gam[:, i]'
    cov_ha[:, i] = vec(tril(tmp))[vec(tril(tmp).!=0.0)]
    big_cov_ha[:, :, i] = Hermitian(tmp, :L)
  end
  (big==true) && return big_cov_ha
  return cov_ha
end
#
function random_gen_delta(t1, t2; seed=89871924)
  srand(seed)
  rnd_delta = rand(Normal(), t1, t2)
  return rnd_delta
end
