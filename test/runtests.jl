using ValueOrientedRiskManagementInsurance
using Base.Test
using DataFrames

t_0 = 0
T = 5
## pricing assumptions
rfr_price = fill(0.005, T)                         ## discount rate for pricing
prob_price = DataFrame( ## biometric probabilities
                       qx = 0.001 .+ [0:(T-1)] * 0.0001,
                       sx = fill(0.1, T)
                       )
λ_price =  DataFrame(  ## cost profile
                     boy = [0.05, rep(0.0, T-1)], ## costs beginning of year
                     eoy = fill(0.04, T),         ## costs end of year
                     infl = fill(0.01, T)         ## cost inflation per year
                     )
β = DataFrame( ## profile insurance contract
              qx = fill(1.0, T),                  ## death benefit rel to insured sum
              sx = cumsum(fill(0.9, T)),          ## lapse ben rel to premium
              px = [fill(0.0, T-1), 1],           ## life ben  rel to insured sum
              prem = fill(1.0, T)                 ## normalized premium payments
              )

product = Product(rfr_price, prob_price, β, λ_price)


## Test of the pricing formula
tst, tst_disc, tst_infl, tst_lx = 0.0, 1.0, 1.0, 1.0
for t in 1:product.dur
  if t > 1
    tst_lx *= product.prob[t-1, :px]
  end
  tst += tst_disc * tst_lx * product.prem_norm
  tst -= tst_disc * tst_lx * product.λ[t, :boy]
  tst_disc /= (1 + product.rfr[t])
  tst_infl *= (1 + product.λ[t, :infl])
  tst -=
    tst_disc * tst_lx * product.prob[t,:qx] * product.β[t, :qx]
  tst -=
    tst_disc * tst_lx * product.prob[t,:px] * product.β[t, :px]
  tst -=
    tst_disc * tst_lx * product.prob[t,:sx] * product.β[t, :sx] *
    product.prem_norm
  tst -= tst_disc * tst_lx * tst_infl * product.λ[t,:eoy]
 end
@test_approx_eq_eps(tst, 0.0, 1e-15)

println("test completed")
