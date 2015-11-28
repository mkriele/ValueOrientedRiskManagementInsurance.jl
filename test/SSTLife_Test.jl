using ValueOrientedRiskManagementInsurance

using Distributions
using DataFrames
include("SSTLife_Input.jl")
include("SSTLife.jl")

println("Start SSTLife test")

@test_approx_eq(round(rtk_start, 2), 158.58)
@test_approx_eq(round(tc, 2), 147.60)
@test_approx_eq(round(sst_ratio, 4), 1.0744)

println("End SSTLife test")
