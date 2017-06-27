# auxiliary simulation functions
#
@doc """
Simulates a stochastic process {yₜ} following an ARMA(1,1) with parameters φ and Θ for t = 1, ..., T periods.

```math
  y_t = \\mu_0 + \\mu_1 + \\phi*y_{t-1} + \\epsilon_t + \\theta*\\epsilon_{t-1}, \\quad \\epsilon_t \\sim \\mathbf{N}(0, \\sigma^2)
```

Required Input Parameters
----------------------------------
none

Optional Input Parameters
----------------------------------
`mu::Float64`: AR(1) process constant term. (default=0.0)

`lin_trend::Float64`: AR(1) process linear time trend. (default=0.0)

`ar::Float64`: AR(1) autorregressive parameter. (default=0.5)

`ma::Float64`: MA(1) autorregressive parameter. (default=0.5)

`T::Int64`: Length of simulated time series. (default=100)

`y_0::Float64`: Initial value of the series. (default=0.0)

`stderr::Float64` : Normally distributed innovation standard error. (default=1.0)

Output
----------------------------------
`y::Vector{Float64}`: Simulated ARMA(1,1) time series of length T.

""" ->
function simul_arma(;mu::Float64 = 0., lin_trend::Float64 = 0., ar::Float64 = .5, ma::Float64 = .5, T::Integer = 100, y_0::Float64 = 0.0, stderr::Float64 = 1.0)
  #y_t = μ + ϕy_t-1 + u_t/ u_t = ɛ_t + Θɛ_t-1
  y = zeros(T); ɛ = rand(Normal(0.0,stderr),T); u = zeros(T);
  for t in 2:T
    u[t, :] = ɛ[t, :] + ma*ɛ[t-1, :] #MA
    y[t, :] = mu + lin_trend*t + ar*y[t-1, :] + u[t, :] #AR
  end
  return y
end
#
#
function simul_llm(; T=100::Int64, x0=0.0::Float64)
  # x_t = μ_t + ϵ_t, ϵ_t∼N(0, σ2_ϵ)
  # μ_t = μ_{t-1} + η_t, ϵ_t∼N(0, σ2_η)
  #state and signal initalized at zero
  x = zeros(T)
  μ = zeros(T)
  σ2_ϵ = 1.0;
  σ2_η = 1.0;
  ϵ = rand(Normal(0.0, σ2_ϵ), T)
  η = rand(Normal(0.0, σ2_η), T)
  x[1] = x0
  for t in 2:T
    μ[t] = μ[t-1] + η[t]
    x[t] = μ[t] + ϵ[t]
  end
  return x
end
