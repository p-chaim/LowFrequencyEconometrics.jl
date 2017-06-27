function mw_h_cov_ha_i0_i1_i2(psi, param, nsim, rmatrix_strig)
#mw_h_cov_ha_i0_i1_i2 - One line description of what the function or script performs (H1 line)
# Generate nsim HA covariance mtrices for h-test for
#   (1). I(0) model.
#   (2). RW model.
#   (3). Integrated RW Model
#   Three covariance matrices are returned.
#   These are the ingredients to the covariance matrices for the local level and integrated local level models.
#
# Syntax:  [output1,output2] = function_name(input1,input2,input3)
#
# Inputs:
#   psi::Array{Float64,2} - dct transform.
#   param::Float64 - local parameters.
#   nsim::Int64 - number of simulations.
#   rmatrix_string = standard random variables.
#
# Outputs:
#   cov_ha_i0 - MATRIX version of cov_ha.
#   cov_ha_i1 - MATRIX version of cov_ha for summed process.
#   cov_ha_i2 - MATRIX version of cov_ha for double summed process.
#
#
#
#
  t = size(psi, 1)
  q = size(psi, 2)
  #qq = Int(q*(q+1)/2)
  # Get draws of variance function
  w = rmatrix_strig
  w = cumsum(w[1:t, 1:nsim], 1)/sqrt(t) # Random Walk
  w = w*(param/sqrt(q)) # Scaled
  sdpath_mat = exp.(w)
  # Construct summation matrix
  sum_mat = zeros(t, t)
  for i ∈ 1:t
    sum_mat[i, 1:i] = ones(1, i)
  end
  sum_mat2 = sum_mat*sum_mat # Double integration
  cov_ha_i0 = NaN*zeros(q, q, nsim)
  cov_ha_i1 = cov_ha_i0
  cov_ha_i2 = cov_ha_i1
  psi_0 = psi'
  psi_1 = psi'sum_mat
  psi_2 = psi'sum_mat2
  for i ∈ 1:nsim
    tmp = psi_0.*sdpath_mat[:,i]'
    tmp = tmp*tmp'
    cov_ha_i0[:, :, i] = Hermitian(tmp, :L)
    tmp = psi_1.*sdpath_mat[:,i]'
    tmp = tmp*tmp'
    cov_ha_i1[:, :, i] = Hermitian(tmp, :L)
    tmp = psi_2.*sdpath_mat[:,i]'
    tmp = tmp*tmp'
    cov_ha_i2[:, :, i] = Hermitian(tmp, :L)
  end
return cov_ha_i0, cov_ha_i1, cov_ha_i2
end
