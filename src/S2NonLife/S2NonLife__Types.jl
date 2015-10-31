export NLLob, PY, CY, NY

const PY, CY, NY = 1, 2, 3   ## previous, current, next year

"""
Line of Business (Non-Life)

  - `name::String`: Name of the line of business
  - `index::Int`: Identifying index of the line of business
  - `prem_gross_w::Vector{Float64}`: Written gross premium
  - `prem_w::Vector{Float64}`: Written premium
  - `prem_gross_cy::Float64`: Gross premium earned
  - `prem::Vector{Float64}`: Premium earned
  - `upr_gross::Vector{Float64}`: Unearned gross prem. reserve
  - `upr::Vector{Float64}`: Unearned premium reserve
  - `re_prop_q::Vector{Float64}`: Prop. reinsurence ceded
  - `β::Vector{Float64}`: Reserve pattern
  - `vol_prem::Float64`:  Volume factor premium
  - `vol_prem_vc::Float64`: Var. coeff. net prem risk
  - `vol_res::Float64`: Volume factor reserve
  - `vol_prem_res_sd::Float64`: Standard deviation prem+res
"""
type NLLob
  name::String
  index::Int
  prem_gross_w::Vector{Float64}
  prem_w::Vector{Float64}
  prem_gross_cy::Float64
  prem::Vector{Float64}
  upr_gross::Vector{Float64}
  upr::Vector{Float64}
  re_prop_q::Vector{Float64}
  β::Vector{Float64}
  vol_prem::Float64
  vol_prem_vc::Float64
  vol_res::Float64
  vol_prem_res_sd::Float64
end

