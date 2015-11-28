using ValueOrientedRiskManagementInsurance

using Distributions
using DataFrames
using Base.Test

include("S2NonLife.jl")

println("Start S2NonLife test")

@test_approx_eq(round(scr_prem_res, 2), 210.95)
@test_approx_eq(round(scr_lapse, 2), 0.0)
@test_approx_eq(round(scr_cat_fire, 2), 60.00)
@test_approx_eq(round(scr_cat_liab, 2), 239.65)
@test_approx_eq(round(scr_cat_man_made, 2), 247.05)
@test_approx_eq(round(scr_cat_other, 2), 0.0)
@test_approx_eq(round(scr_cat_nprop, 2), 0.0)
@test_approx_eq(round(scr_cat, 2), 247.05)
@test_approx_eq(round(scr_nl, 2), 362.76)

println("Start S2NonLife test")
