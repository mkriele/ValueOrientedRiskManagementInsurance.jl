using ValueOrientedRiskManagementInsurance
using DataFrames
using Base.Test

include("ECModel.jl")

println("Start ECModel test")

@test_approx_eq(round(total.gross.profit_mean, 2), 293.08)
@test_approx_eq(round(total.net.profit_mean, 2), 249.39)
@test_approx_eq(round(total.gross.eco_cap, 2), 1096.71)
@test_approx_eq(round(total.net.eco_cap, 2), 823.70)
@test_approx_eq(round(total.gross.rorac, 3), 0.267)
@test_approx_eq(round(total.net.rorac, 3), 0.303)
## risk adjusted pricing:
@test_approx_eq(round(ins_input_rp[1, :loss_ratio], 4), 0.7684)
@test_approx_eq(round(ins_input_rp[2, :loss_ratio], 4), 0.7158)
@test_approx_eq(round(ins_input_rp[3, :loss_ratio], 4), 0.7500)
## capital optimization
@test_approx_eq(round(rorac_net_oc[i_oc_opt], 4), 0.2690)
## product mix optimization, chosen quote f
@test_approx_eq(round(avg_ceded_ofr, 4), 0.4928)
@test_approx_eq(round(rorac_net_or_opt, 4), 0.3376)

println("End ECModel test")
