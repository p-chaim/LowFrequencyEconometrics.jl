function ar1_cov(T, rho)
  acv = rho.^(0:T-1);
  acv = acv/(1-rho^2);
  covmat = acv2cov(acv)
  corr = acv./acv[1]
  return covmat, corr
end
