function den_invariant(x,sigx_inv);
  # Compute Density --- scale included
  m = size(x,1);
  xc = chol(Hermitian(sigx_inv))*x;
  #xc = factorize(sigx_inv).factors*x;
  c = sqrt(det(sigx_inv))*0.5*gamma(m/2)/(pi^(m/2));
  den = c*((sum(xc.^2,1)).^(-m/2))';
  return den
end
