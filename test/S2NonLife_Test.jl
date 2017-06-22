using ValueOrientedRiskManagementInsurance

using Distributions
using DataFrames
using Base.Test

include("S2NonLife.jl")

println("Start S2NonLife test")

@test round(scr_prem_res, 2) ≈ 210.95
@test round(scr_lapse, 2) ≈ 0.0
@test round(scr_cat_fire, 2) ≈ 60.00
@test round(scr_cat_liab, 2) ≈ 239.65
@test round(scr_cat_man_made, 2) ≈ 247.05
@test round(scr_cat_other, 2) ≈ 0.0
@test round(scr_cat_nprop, 2) ≈ 0.0
@test round(scr_cat, 2) ≈ 247.05
@test round(scr_nl, 2) ≈ 362.76

println("End S2NonLife test")
