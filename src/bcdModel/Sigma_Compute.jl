#
function q_Sigma_cd_grid(q::Int, c_grid::Array{Float64,1}, d_grid::Array{Float64,1})
    mkΣ_T = Int(1000)#numerical integration grid
    Tmax = mkΣ_T + 2
# Compute cosine transforms and weights
  fvec = collect(π*(1:1:q))'
  tvec = collect((0.5:1:mkΣ_T)/mkΣ_T)
  v = zeros(Tmax, q)
  v[1:mkΣ_T, 1:q] = sqrt(2)*cos(tvec*fvec)
  v_1 = v
#
  x = (1/mkΣ_T)*collect(1:Tmax-1)
  n_c = length(c_grid)
  n_d = length(d_grid)
  qΣcd_grid = NaN*zeros(n_c, n_d, q, q)
    #
  for ic in 1:n_c
    println("c iter $ic of $n_c")
    c = c_grid[ic]
    (c == 0) && (c = 0.01) # corner case
    cx = c*x
    # Summation by parts formula used in d0 < 0.5;
    vqd = NaN*zeros(Tmax, q)
    vqd[1:Tmax-1, :] = (c/mkΣ_T)*v_1[1:Tmax-1, :] - (v_1[2:Tmax, :] - v_1[1:Tmax-1, :])
    vqd[1, :] = vqd[1, :] - v_1[1, :]
    vqd[Tmax, :] = 0
    v_0 = mkΣ_T*vqd
    for id in 1:n_d
      d = d_grid[id]
      vqd = NaN*zeros(Tmax, q)
      d0 = d;
      (d < 0.5) && (d = d+1)
      (d == 0.5) && (d = 0.51)
      if (d0 < 0.5)
        v = v_0
      else
        v = v_1
      end
      #
      ga_0 = c^(1-2*d)*sqrt(π)*gamma(-0.5+d)/gamma(d)
      bselval = besselk(d-0.5, cx)
      ga = ((2^(1.5-d))*sqrt(π)/gamma(d))*((x/c).^(d-0.5)).*bselval
      vqd[1, :] = [ga_0 ga']*v
      for t in 1:Tmax-2
        vqd[t+1, :] = [flipdim(ga[1:t], 1)' ga_0 ga[1:Tmax-t-1]']*v
      end
      vqd[Tmax, :] = [flipdim(ga, 1)' ga_0]*v
      qΣcd = v'*vqd/(mkΣ_T^2)
      qΣcd = qΣcd/(2*π)
      #Normalize to have average value equal to unity
      tmp = mean(diag(qΣcd))
      qΣcd = qΣcd/tmp
      qΣcd_grid[ic, id, :, :] = qΣcd
    end
  end
  return qΣcd_grid
end
#
function interp_q_Sigma_cd_grid(qΣcd_grid, c, d, c_grid, d_grid)
  q = size(qΣcd_grid, 4)
  qΣcd = NaN*zeros(q,q)
  for i in 1:q
    for j in 1:q
      tmp = qΣcd_grid[:, :, i, j] # these are the values of the covariance between transforms i and j for each value of c and d
      spl = Spline2D(c_grid, d_grid, tmp)
      qΣcd[i, j] = evalgrid(spl, [c], [d])[1]
    end
  end
  tmp = mean(diag(qΣcd))
  qΣcd = qΣcd/tmp
  return qΣcd
end
#
function estimate_d(lr::bcdModel)
  q = lr.q
  cvec = lr.c.grid
  dvec = lr.d.grid
  Sigma_d_mat=NaN*zeros(lr.q, lr.q, length(dvec))
  Sigma_inv_d_mat = Sigma_d_mat
  #"Estimate" d via grid through likelihood calibration
  lvec = NaN*zeros(length(dvec))
  X_norm = lr.cost
  for i in 1:length(dvec)
    tic()
    d = dvec[i]
    Sigma = Sigma_Compute2(0, 0, d, lr.q)
    sig_inv = inv(Sigma)
    lvec[i, :] = den_invariant(X_norm, sig_inv)
    Sigma_d_mat[:, :, i] = Sigma
    Sigma_inv_d_mat[:, :, i] = sig_inv
    toc()
  end
  lr_stat = mean(lvec)./lvec
  pct_vec = 1-[2/3 0.9 0.95]'
  lr_d_cv = NaN*zeros(length(dvec), length(pct_vec))
  nsim = 50000
  Ysim = randn(lr.q, nsim)
  lvec_sim = NaN*zeros(size(dvec,1),nsim);
  for i in 1:length(dvec)
    tic()
    sig_ho = Sigma_d_mat[:, :, i]
    chol_ho = chol(Hermitian(sig_ho))
    Xsim = chol_ho*Ysim
    tmp = sqrt(sum(Xsim.^2, 1))
    Xsim_norm = Xsim./repmat(tmp,q,1);
    lvec_sim = NaN*zeros(size(dvec,1),nsim);
    for j = 1:size(dvec,1);
         sig_inv = Sigma_inv_d_mat[:,:,j];
         lvec_sim[j,:] = den_invariant(Xsim_norm,sig_inv);
     end;
     lr_stat_ho = mean(lvec_sim)./lvec_sim[i,:];
     tmp = pctile(lr_stat_ho',pct_vec);
     lr_d_cv[i, :] = tmp';
     toc()
   end
   tmp, ii = findmax(lvec)
   println("MLE (from grid) for d: $(dvec[ii])")
   condmat = repmat(log(lr_stat), 1, size(lr_d_cv, 2)) .> log(lr_d_cv)
   println("Confidence sets for d")
   cset = repmat(pct_vec, 1, 2)
   for i in 1:length(pct_vec)
     tmp = dvec[condmat[:, i].==true]
     #println(" Level = $(pct_vec[i]): $(tmp[1]) $(tmp[end])")
     #cset[i, :] = [tmp[1] tmp[end]]
   end
   #tmp, ii = findmin(abs(dvec-0))
   #lvec = lvec
   lr.d.val = dvec[ii]
   lr.stuf[:lvec] = lvec
   lr.d.lik = lvec
  #return lvec, lr_stat, Sigma_d_mat, Sigma_inv_d_mat, lr_d_cv, lvec_sim, cset
  return lr
end
#
function Sigma_Compute2(b, c, d, q; iconst = false)
  T_general = 150; T_special = 600;
  cov_Id = 0
  if c < 0.0001
    T = T_special;
    cov_Id = frac_cov(T,d)
  else
    if d == 1
      T = T_special;
      rho = 1 - c/T;
      rho = max(0, rho);
      cov_Id, tmp = ar1_cov(T,rho);
    else
      error("not implemented")
      T = T_general;
      rho = 1 - c/T;
      rho = max(0, rho);
      #cov_Id = frac_dif_cov(T, rho, d);
      cov_Id, tmp = ar1_cov(T,rho);
    end
  end
  #Add I(0) components
  scale_I0 = (b*(T^d))^2;
  #println(cov_Id); println(scale_I0); println(eye(T))
  covmat = cov_Id + scale_I0*eye(T);
  psi = psi_compute(T, q);
  if iconst == true
    psi = [ones(T,1)/T psi];
  end
  # Step 3: compute Sigma
  Sigma = psi'*covmat*psi;
  expoente = 1. - 2*d;
  Sigma = (T^expoente)*Sigma;
  #= Verify Sigma is positive definite
  ier = 0;
  tmp, tmp2 = eig(Sigma);
  if minimum(tmp) < 0
    ier = 1;
    Sigma_2 = Sigma*Sigma';
    V, D2 = eig(Sigma_2);
    D_vec = sqrt(diag(D2));
    D_mat = diagm(D_vec);
    Sigma = V*D_mat*V';
  end =#
  return Sigma
end
#
