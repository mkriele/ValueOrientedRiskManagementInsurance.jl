using ValueOrientedRiskManagementInsurance
using DataFrames
using Base.Test

include("ECModel.jl")

println("Start ECModel test")

@test round(total.gross.profit_mean, 2) ≈ 293.08
@test round(total.net.profit_mean, 2) ≈ 249.39
@test round(total.gross.eco_cap, 2) ≈ 1096.71
@test round(total.net.eco_cap, 2) ≈ 823.70
@test round(total.gross.rorac, 3) ≈ 0.267
@test round(total.net.rorac, 3) ≈ 0.303
## risk adjusted pricing:
@test round(ins_input_rp[1, :loss_ratio], 4) ≈ 0.7684
@test round(ins_input_rp[2, :loss_ratio], 4) ≈ 0.7158
@test round(ins_input_rp[3, :loss_ratio], 4) ≈ 0.7500
## capital optimization
@test round(rorac_net_oc[i_oc_opt], 4) ≈ 0.2690
## product mix optimization, chosen quote f
@test round(avg_ceded_ofr, 4) ≈ 0.4928
@test round(rorac_net_or_opt, 4) ≈ 0.3376

println("End ECModel test")
