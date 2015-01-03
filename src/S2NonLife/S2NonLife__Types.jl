export NLLob, PY, CY, NY

const PY, CY, NY = 1, 2, 3   ## previous, current, next year

#################################################################

type NLLob  ## Non-life line of business
  name::String
  index::Int
  prem_gross_w::Vector{Float64}  ## written gross premium
  prem_w::Vector{Float64}        ## written premium
  prem_gross_cy::Float64         ## gross premium earned
  prem::Vector{Float64}          ## premium earned
  upr_gross::Vector{Float64}     ## unearned gross prem. reserve
  upr::Vector{Float64}           ## unearned premium reserve
  re_prop_q::Vector{Float64}     ## prop. reinsurence ceded
  β::Vector{Float64}             ## reserve pattern
  vol_prem::Float64              ## volume factor premium
  vol_prem_vc::Float64           ## var. coeff. net prem risk
  vol_res::Float64               ## volume factor reserve
  vol_prem_res_sd::Float64       ## standard deviation prem+res
end

