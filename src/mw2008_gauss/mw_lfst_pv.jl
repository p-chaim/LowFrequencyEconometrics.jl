function  mw_lfst_pv(x_t, itrend, param_m, param_t, nsim)
# mwlfst_pv - rewrite of Muller and Watson (2015) GAUSS code. Calculates LFST
#   test pvalues.
# Syntax:  pv::Float64 = mwlfst_pv(x_t, itrend, param_m, param_t, nsim)
#
# Inputs:
#   x_t::Vector{Float64} - Discrete cosine transforms of data.
#   itrend::Int64 - 1 if x_t is demeaned, 2 if x_t is detrended.
#   param_m::Float64 - "g" parameter for demeaned alternative.
#   param_t::Float64 - "g" parameter for detrended alternative.
#   nsim::Int64 - Number of simulations to compute p-value
# Outputs:
#   pv::Float64 - p-value
#------------- BEGIN CODE --------------
t = 1000 # Sample Size Used for Asymptotic Approximation
q = size(x, 1)
# Compute Null and Alternative Covariance Matrices
if itrend == 1
  psi = psi_compute(t, q)
  cov_h0 = psi'psi
  g_t = param_m/t
  gt2 = g_t.^2
  tmp = mwrw_cov(t)
  tmp = gt2*tmp
  tmp1 = diag(tmp) + ones(t, 1)
  cov_u = diagrv(tmp, tmp1)
  cov_ha = psi'cov_u*psi
end
# Compute Test Statistic @
  chol_cov_h0 = chol(cov_h0)
  a = inv(chol_cov_h0')
  x_t_trans = a*x_t # @ X transformed under Ho @
  cov_ha_trans = a*cov_ha*a' # @ Covariance under Ha Transformed @
#   @ Note ... LF under Ho is unity @
  lr = mw_marg_lf(x_t_trans, cov_ha_trans)
# @ Simulations for Pvalue @
  xsim = rand(Normal(), size(x,1), nsim)
  lr_sim = mw_marg_lf(xsim, cov_ha_trans)
  pv = mean(lr_sim .>(lr - 0.0001))
return lr, pv
end

function diagrv(x, v)
  y = x - diagm(diag(x)) + diagm(v[:])
  return y
end
