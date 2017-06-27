function mw_rw_cov(t)
  sum_mat = zeros(t,t)
  for i ∈ 1:t
    sum_mat[i, 1:i] = ones(1, i)
  end
  cov = sum_mat*sum_mat'
  return cov
end
