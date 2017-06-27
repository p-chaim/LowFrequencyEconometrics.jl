module LowFrequencyEconometrics

#Import your packages
using KernelDensity,
		Distributions,
		Dierckx,
		StatsBase,
		DataFrames,
		PyPlot
#import Base: func1 #Any function you add dispatches to need to be imported directly

#abstract AbType #Define abstract types before the types they abstract!
include("util\\psi_compute.jl")
include("util\\simul_.jl")
include("util\\den_invariant.jl")
include("util\\frac_cov.jl")
include("util\\acv2cov.jl")
include("util\\ar1_cov.jl")
include("bcdModel\\bcdModel.jl") #Include all the functionality
include("bcdModel\\bcd_spectrum.jl") #Include all the functionality
include("bcdModel\\Sigma_Compute.jl") #Include all the functionality
include("tests\\lfH_test.jl")
include("tests\\lfS_test.jl")
include("tests\\lfst_test.jl")
include("tests\\lfur_test.jl")

export 
#from "util\\simul_.jl"
simul_arma,
simul_llm,
#from "util\\frac_cov.jl"
frac_cov,
#from "util\\psi_compute.jl"
psi_compute,
#from "util\\den_invariant.jl"
den_invariant,
#from "util\\acv2cov.jl"
acv2cov,
#from "util\\ar1_cov.jl"
ar1_cov,
#from "bcdModel\\bcd_spectrum.jl"
plot_locS, #Export the functions you want users to use
#from "bcdModel\\bcdModel.jl"
bcdModel,
#from "bcdModel\\Sigma_Compute.jl"
q_Sigma_cd_grid,
interp_q_Sigma_cd_grid,
estimate_d,
Sigma_Compute2,
#from "tests\\lfH_test.jl"
#from "tests\\lfS_test.jl"
#from "tests\\lfst_test.jl"
lfst_test,
lfst_pv,
mw_lfst_pv
#from "tests\\lfur_test.jl"

end
