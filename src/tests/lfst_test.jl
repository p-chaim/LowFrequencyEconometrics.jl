@doc """
```math
LFST = \\frac{\\sum^q_{j=1} X^2_{j,T}}{\\sum^q_{j=1}}
```
""" ->
function lfst_test(X, g)
  q = size(X, 1)
  tmp = (collect(1:1:q)*π).^2
  cov_rw = 1./tmp
  cov_ha = ones(q).+(g*g*cov_rw)
  cov_ha_inv = 1./cov_ha
  chol_ha = sqrt(cov_ha_inv)
  Xa = X.*chol_ha
  lfst = sum(X.^2)/sum(Xa.^2)
  nsim = 50000
  Xsim = randn(q, nsim)
  Xasim = Xsim.*repmat(chol_ha, 1, nsim)
  nx = sum(Xsim.^2, 1)
  dx = sum(Xasim.^2, 1)
  lfst_sim = (nx./dx)'
  lfst_pvalue = mean(lfst_sim.>=lfst)
  return lfst, lfst_pvalue
end
function lfst_pv(X, g)
  q = size(X, 1)
  tmp = (collect(1:1:q)*π).^2
  cov_rw = 1./tmp
  cov_ha = ones(q).+(g*g*cov_rw)
  cov_ha_inv = 1./cov_ha
  chol_ha = sqrt(cov_ha_inv)
  Xa = X.*chol_ha
  lfst = sum(X.^2)/sum(Xa.^2)
  nsim = 50000
  Xsim = randn(q, nsim)
  Xasim = Xsim.*repmat(chol_ha, 1, nsim)
  nx = sum(Xsim.^2, 1)
  dx = sum(Xasim.^2, 1)
  lfst_sim = (nx./dx)'
  lfst_pvalue = mean(lfst_sim.>=lfst)
  return lfst, lfst_sim
end

function  mw_lfst_pv(x_t, itrend, param_m, param_t, nsim; stat=false)
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
q = size(x_t, 1)
# Compute Null and Alternative Covariance Matrices
if itrend == 1
  psi = psi_compute(t, q)
  cov_h0 = psi'psi
  g_t = param_m/t
  gt2 = g_t.^2
  tmp = mw_rw_cov(t)
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
  xsim = rand(Normal(), size(x_t,1), nsim)
  lr_sim = mw_marg_lf(xsim, cov_ha_trans)
  pv = mean(lr_sim .>(lr - 0.0001))
return lr, lr_sim
end

function diagrv(x, v)
  y = x - diagm(diag(x)) + diagm(v[:])
  return y
end

function mw_marg_lf(x::Vector{Float64}, cov::Array{Float64,2})
# This is a scaled version of the density in King and Hillier (1985) Equation (2) @
# -- Note x is ncov x n ... compute n marginal likelihood values @
# Compute Marginal Likelihood value for x with covariance matrix cov
  x = x./sqrt(sum(x.^2, 1))'
  m = size(x, 1)
  d = det(cov)
  # (gauss code) x = chol(invpd(cov))
  c = chol(Hermitian(inv(x'x)))
  #x = c*x
  tmp = sum(x, 1)
  tmp = d*(tmp.^m)
  mlf = sqrt(1.0./tmp)
  return mlf
end
function mw_marg_lf(x::Array{Float64,2}, cov::Array{Float64,2})
  nsim = size(x, 2)
  mlf_vec = Vector{Float64}(nsim)
  for isim ∈ 1:nsim
    mlf_vec[isim, :] = mw_marg_lf(x[:, isim], cov)
  end
return mlf_vec
end

function mw_rw_cov(t)
  sum_mat = zeros(t,t)
  for i ∈ 1:t
    sum_mat[i, 1:i] = ones(1, i)
  end
  cov = sum_mat*sum_mat'
  return cov
end
