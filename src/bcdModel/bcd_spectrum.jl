
function locS(w::Real, b::Real, c::Real, d::Real)
    S = (w^2 + c^2)^(-d) + b^2
  return S
end

function locS(b::Real, c::Real, d::Real)
  S = NaN*zeros(10000); n_w = length(S)
  w_grid = collect(linspace(0.001, 15*π, n_w))
  i = 0
  for w in w_grid
    i += 1
    S[i] = (w^2 + c^2)^(-d) + b^2
  end
  return S
end


function locSs(b::Real, c::Real, d::Real)
  S = locS(b, c, d)
  n = length(S)
  w_grid = collect(linspace(0.001, 15π, n))
  tmp, i1 = findmin(abs(w_grid-6*π))
  Ss = S./S[i1]
  x = w_grid/π
  return Ss, x
end
#
function plot_locS(b::Real, c::Real, d::Real)
    title("Local to Zero Spectrum: S(ω;b,c,d)∝(c² + d²)^(-1) + b²")
    Ss, x = locSs(b, c, d)
    _plt = plot(x, log(Ss), label = "($(round(b, 2)), $(round(c, 2)), $(round(d, 2)))")
    legend()
    return _plt
end
#
function plot_locS(m::String)
  if m == "local level"
    println("$m model")
    b_grid = [1/50 1/10 1/2 10e16]
    for b in b_grid
      Ss, x = locSs(b, 0, 1)
      plot(x, log(Ss), label = "b = $b")
      legend()
    end

  end
  if m == "local to unity"
    println("$m model")
    c_grid = [0 3 10 30]
    for c in c_grid
      Ss, x = locSs(0, c, 1)
      plot(x, log(Ss), label = "c = $c")
    end
    legend()
  end
  if m == "fractional"
    println("$m model")
    d_grid = [-1/3 0 1/3 1 4/3]
    for d in d_grid
      Ss, x = locSs(0, 0, d)
      plot(x, log(Ss), label = "d = $d")
    end
    legend()
  end
  if m == "bcd"
    b_grid = [1/50 1/10 1/2 10e16]
    c_grid = [0 3 10 30]
    d_grid = [-1/3 0 1/3 1 4/3]
    for b in b_grid
      for c in c_grid
        for d in d_grid
          println("$b, $c, $d")
          Ss, x = locSs(b, c, d)
          plot(x, log(Ss))
        end
      end
    end
  end
end

#end #module
