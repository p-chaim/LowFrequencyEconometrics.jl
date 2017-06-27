function mw_marg_lf(x::Vector{Float64}, cov::Array{Float64,2})
# This is a scaled version of the density in King and Hillier (1985) Equation (2) @
# -- Note x is ncov x n ... compute n marginal likelihood values @
# Compute Marginal Likelihood value for x with covariance matrix cov
  x = x./sqrt(sum(x.^2, 1))'
  m = size(x, 1)
  d = det(cov)
  # (gauss code) x = chol(invpd(cov))
  c = chol(Hermitian(inv(x'x)))
  print(x)
  #x = c*x
  tmp = sum(x, 1)
  tmp = d*(tmp.^m)
  mlf = sqrt(1.0./tmp)[:]
  return mlf
end
function mw_marg_lf(x::Array{Float64,2}, cov::Array{Float64,2})
  nsim = size(x, 2)
  println(nsim)
  mlf_vec = Vector{Float64}(nsim)
  for isim ∈ 1:nsim
    mlf_vec[isim] = mw_marg_lf(x[:, isim], cov)[1]
  end
return mlf_vec
end
