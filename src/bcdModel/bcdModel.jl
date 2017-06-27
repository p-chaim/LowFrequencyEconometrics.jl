# Univariate bcdModel type and constructors.
type bcdPar
  val::Float64
  grid::Vector{Float64}
  lik::Vector{Float64}
  lik_max_index::Int64
  function bcdPar(val::Float64, grid::Vector{Float64})
    self = new()
    self.val = val
    self.grid = grid
    lik = NaN*zeros(length(self.grid))
    lik_max_index = 0
    return self
  end
end


type bcdModel
  # Class Attribute
    name::String # Series name
    vals::Vector{Float64}#Raw series
    dte::Vector{Date} # series dates in Date format
    q::Int64 # number of cosine transforms
    T::Int64
    avg::Float64
    lrstd::Float64
    std::Float64
    cost::Vector{Float64,} # low frequency cosine transforms
    scost::Array{Float64,} # standardized low frequency cosine transforms
    proj::Vector{Float64} #
    Omega::Array{Float64,} # Sample covariance matrix in (qxq) form
    #parameterization
    b::bcdPar
    c::bcdPar
    d::bcdPar
    Sigma::Array{Float64,} # qΣbcd model-specific covariance matrix. Note this depends only on parameters b, c, d, and on cosine transform weights cost.
    stuf # dict to dump other shit
  # Class methods
    getCov::Function
    mle_d::Function
  # bcdModel Constructor
  function bcdModel(vals::Array{Float64,1}; q=12::Int64, b=0.0::Float64, c=5.0::Float64, d=0.3::Float64, name="Series"::String, start_dte=1900, nc = 5, nd = 50, cov_calc=false)
    #####################
    # Class constructor #
    #####################
    self = new()
    self.name = name
    self.vals = vals
    self.dte = collect(Date(start_dte):Dates.Year(1):Date(start_dte+length(vals)-1))
    self.q = q
    self.T = length(self.vals)
    self.avg = mean(self.vals)
    self.std = std(self.vals)
    #self.cost = NaN*zeros(q)
      Ψ = psi_compute(self.T, self.q)
    self.cost = Ψ'*(self.vals - self.avg)
    self.lrstd = sqrt(mean(self.cost.^2)*length(self.vals))
    self.scost = (self.cost)./sqrt(self.cost'*self.cost)
    #self.scost = NaN*similar(self.cost)
    #self.proj = NaN*similar(self.vals)
    self.proj = self.avg + (Ψ/(Ψ'Ψ))*self.cost
    self.Omega = NaN*zeros(q, q)
    self.b = bcdPar(b, [0.])
      (nc <= 1) && (nc = nc + 1)
      c_grid = zeros(nc)
      for ic in 2:nc
        c_grid[ic] = 2.*(200)^(2.^(ic./(nc .- 1.0)))
      end
    self.c = bcdPar(c, c_grid)
      (nd <= 1) && (nd = nd + 1)
      d_grid = collect(linspace(-0.49, 1.49, nd))
    self.d = bcdPar(d, d_grid)
    self.stuf = Dict()
    #################
    # Class methods #
    #################
    self.getCov = function()
      Sigma_grid = q_Sigma_cd_grid(self.q, self.c.grid, self.d.grid)
      Sigma = interp_q_Sigma_cd_grid(Sigma_grid, self.c.val, self.d.val, self.c.grid, self.d.grid)
      self.Sigma = Sigma
      return self
    end
    #
    self.mle_d = function()
      Sigma_d_mat = NaN*zeros(self.q, self.q, length(self.d.grid))
      dcount = 0.
      self.d.lik = NaN*similar(self.d.grid)
      for dlik in 1:length(self.d.grid)
        tic
        Sigma = Sigma_Compute2(0.0, 0.0, self.d.grid[dlik], self.q)
        Sigma_inv = inv(Sigma)
        self.d.lik[dlik] = den_invariant(self.scost, Sigma_inv)[1]
#        Sigma_d_mat[:, :, dlik] = Sigma
#        Sigma_inv_d_mat[:, :, dlik] = Sigma_inv
        tic()
      end
      tmp, dmax_index = findmax(self.d.lik)
      self.d.val = self.d.grid[dmax_index]
      self.d.lik_max_index = dmax_index
      println("d_MLE: $(self.d.val) (estimated from an equally-spaced grid of $(length(self.d.grid)) points between -0.49 and 1.49 for the fractional difference model, i.e. bcd(0, 0, 1))")
      return self
    end
    return self
  end
end
