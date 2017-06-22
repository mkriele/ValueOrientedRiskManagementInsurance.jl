using Distributions
using DataFrames
using Base.Test
using Dierckx
using ValueOrientedRiskManagementInsurance


import ValueOrientedRiskManagementInsurance.zb

### TESTS
curr_dir = dirname(@__FILE__())

println("start HWBS.jl...")

## All Data sourced from Swiss National Bank (www.snb.ch)
# Read the raw data for market prices.  The raw data only
# constains spot rates for different measurement dates
# (columns :xyyyy_mm, where yyyy is the year and mm the month)
# and the corresponding Durations (column :Duration).
df_chf_spots = readtable( curr_dir * "/HWBS_Input_CHF_Spot.csv",
                          header = true)
# contains various short rates ("SARON", "Overnight",
# "MM_3months") for different measurement dates
df_chf_shorts = readtable(curr_dir * "/HWBS_Input_CHF_Short.csv",
                          header = true)
# parameters for the Hull-White Black-Scholes model
hwpar = readtable(curr_dir * "/HWBS_Input_Parameters.csv",
                  header = true)
# Asset portfolio
assets = readtable( curr_dir * "/HWBS_Input_Portfolio.csv",
                    header = true)

# t₀ refers to the time at which we will base the spot curve
# it is given as a column identifyer for df_chf_spots_sparse
t₀ = :x2016_09
# We add the 3month money market debt rate (Switzerland - CHF -
# Money market debt register claims of the Swiss Confederation,
# 3-month) as a proxi for the spot rate of duration 0.
# This is used for extrapolation to spot rates with
# duration < 1 year
spot₀ =
  df_chf_shorts[find(x->x=="SARON",df_chf_shorts[:Type])[1], t₀]

df_prices = DataFrame()
df_prices[:t] = vcat(0.0, df_chf_spots[:Duration])
df_prices[:spot] = vcat(spot₀, df_chf_spots[t₀])

dict_hwpar =
  Dict{Symbol, Real}( Symbol( hwpar[i, :name]) =>
                      hwpar[i, :value] for i = 1:nrow(hwpar))


hw = HWBS(dict_hwpar, df_prices)

asset_id = Array{Symbol}(nrow(assets))
for row in 1:nrow(assets)
  if assets[row,:Type] in ["zb"]
    asset_id[row] =
      assets[row,:Type] *"_" * string(assets[row,:Maturity])
  else
    asset_id[row] = assets[row,:Type]
  end
end

n_scen = 100_000
t₁ = 1
V₀ = zeros(nrow(assets))            ## Value of assets at time 0
V₁ = zeros(nrow(assets), n_scen)    ## Value of assets at time 1
## loss in value from 0 to 1, discounted to time 0
loss = DataFrame()
## random variables for capital market model at time 1
cm₁ = rand_rw(hw, t₁, n_scen)
r₁ = cm₁[1,:]  ## random interest rates at time 1
S₁ = cm₁[2,:]  ## random index values at time 1

discount = zb(hw, 1) ## discount rate from 0 to 1

for row in 1:nrow(assets)
  if assets[row, :Type] == "index"
    V₀[row] = assets[row, :Nominal]
    V₁[row, :] = V₀[row] * S₁
  else
    V₀[row] =
      zb(hw, assets[row, :Maturity]) * assets[row, :Nominal]
    zb_tmp(r) = zb(hw, r, t₁, assets[row, :Maturity])
    V₁[row, :] = zb_tmp.(r₁) * assets[row, :Nominal]
  end
  loss[asset_id[row]] = V₀[row] - discount * V₁[row,:]
end
loss[:total] = sum(loss[asset_id[row]] for row in 1:nrow(assets))

es_assets = es(loss[:total], 0.99)


println("...end HWBS.jl")
