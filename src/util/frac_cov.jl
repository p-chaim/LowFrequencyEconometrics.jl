function frac_cov(t, d)
  dsave = d
  if d < -0.5;
   error("d must be > -0.5");
  end;
  if d > 1.5;
     error("d must be < 1.5");
  end;
  (d>0.50) && (d = d - 1)
  (d == 0) && return covmat = eye(t)
  (abs(d) == 0.5) && (d = 0.49999*sign(d))
  acv = zeros(t, 1)
  acv[1] = 1
  (d > 0) && (Γmd = gamma(d))
  (d < 0) && (Γmd = gamma(1+d)/d)
  Γ_md = gamma(1-d)
  fac = Γ_md/Γmd
  t1 = minimum([100, t-1])
  h = collect(1:t1)
  Γhd = gamma(h+d)
  Γhmd = gamma(h-d+1)
  acv[2:t1+1] = (Γhd./Γhmd)*fac
  if t1+1<t
    hvec = collect(t1+1:t-1)
    anum = lgamma(d+hvec)
    aden = lgamma(hvec-d+1)
    acv[t1+2:t] = fac*exp(anum-aden)
  end
  covmat = acv2cov(acv)
  var_scl=gamma(1-2*d)/(gamma(1-d))^2;
  covmat=covmat*var_scl;
  if dsave > .50
    sum_mat = tril(ones(t,t))
    covmat = sum_mat*covmat*sum_mat'
  end
   return covmat
end
