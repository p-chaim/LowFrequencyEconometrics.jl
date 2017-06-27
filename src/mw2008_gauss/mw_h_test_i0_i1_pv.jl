# mw_htest_i0_i1_pv - Construct H-test pvalue for the I(0) or Unit Root Null Model
# Syntax:  [output1,output2] = function_name(input1,input2,input3)
#
# Inputs:
#   x_t::Vector{Float64} - dct transformed data
#   param::Float64 - parameter for SD in H-test
#   itrend::Int64 - 1 for demean, 2 for detrend
#   ind_int::Int64 - 0 for I(0), 1 for I(1)
#   n_wsim::Int64 - number of simulated sd paths
#   rmatrix_string = name of matrix with standard random variables ... on disk ... w draws
#   nsim = number of simulations for p-value
# Outputs:
#    pv::FLoat64 - pvalue for h test
#
# Example:
#    Line 1 of example
#    Line 2 of example
#    Line 3 of example
#
# Other m-files required: none
# Subfunctions: none
# MAT-files required: none
#
# See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

# Author: Pedro Chaim
# git repository: http://www.
# May 2004; Last revision: 12-May-2004
function mw_h_test_i0_i1_pv(x_t, param, itrend, ind_int, nwsim, rmatrix_string, nsim)
  q = size(x_t, 1)
  qq = q*(q+1)/2
  t = 100 # sample size for asymptotic approximations
  if itrend == 1
    psi = psi_compute(t, q)
  end
  # Construct Appropriate Covariance Matrices @
  cov_ha_i0, cov_ha_i1, cov_ha_i2 = mw_h_cov_ha_i0_i1_i2(psi, param, nwsim, rmatrix_string)
  if ind_int == 0
    cov_h0 = psi'psi
    cov_ha = cov_ha_i0
  elseif ind_int == 1
    cov_u = mw_rw_cov(t)
    cov_h0 = psi'cov_u*psi
    cov_ha = cov_ha_i1
  end
  # Construct Test Statistic @
  lr = mw_lrstat_vec(x_t, cov_h0, cov_ha_i0)
  # @ Compute Simulations for p-value @
  xsim = randn(q, nsim)
  a = chol(Hermitian(cov_h0, :L))
  xsim = a'xsim
  lr_sim = mw_lrstat_vec(xsim, cov_h0, cov_ha)
  pv = mean(lr_sim .> lr, 1)
return pv, lr#, lr, x_t, cov_h0, cov_ha_i0, cov_ha_i1, cov_ha_i2
end
