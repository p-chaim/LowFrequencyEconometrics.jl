function lfur_pv(X, c)
  #
  #
  q = size(X, 1);
  tmp = ((1:1:q)*pi).^2;
  cov_H0 = 1./tmp;
  chol_H0 = sqrt(cov_H0);
  chol_H0_inv = sqrt(1./cov_H0);
  X0 = chol_H0_inv.*X;
  cov_HA = Sigma_Compute2(0, c, 1, q)
  chol_HA_inv = chol(Hermitian(inv(cov_HA), :L));
  XA = chol_HA_inv*X;
  lfur = sum(X0.^2)/sum(XA.^2);
  # Generate 1000 draws under Ho for pvalue calculation
  nsim = 50000;
  Xsim = rand(Normal(0,1), q, nsim).*repmat(chol_H0, 1, nsim);
  X0sim = Xsim.*repmat(chol_H0_inv, 1, nsim);
  XAsim = chol_HA_inv*Xsim;
  nx = sum(X0sim.^2, 1);
  dx = sum(XAsim.^2, 1);
  lfur_sim = squeeze((nx./dx)',2);
  lfur_pvalue =  mean(lfur_sim.>=lfur)
  return lfur, lfur_pvalue
end
