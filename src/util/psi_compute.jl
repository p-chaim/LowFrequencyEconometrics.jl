#
function psi_compute(T::Int64, q::Int64; itrend::Int64=1)
  fvec = π*(1:1:q);
  tvec = (0.5:1:T)/T;
  psi = (sqrt(2))*(cos(tvec*fvec'))./(sqrt(T));
  fvec_2t = fvec./(2*T)
  iota = sin(fvec_2t)./fvec_2t
  if itrend == 1
    psi = psi.*iota'
  end
return psi
end
#
