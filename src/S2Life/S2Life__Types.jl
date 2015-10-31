export GROSS, NET
export S2, S2Mkt, S2MktEq, S2MktInt
export ProjParam

const GROSS, NET = 1,2

abstract S2Module

type S2NotImplemented <: S2Module
  scr::Vector{Float64}
end

function S2NotImplemented(x...)
  return S2NotImplemented(zeros(Float64, 2))
end

typealias S2MktProp S2NotImplemented
typealias S2MktSpread S2NotImplemented
typealias S2MktFx S2NotImplemented
typealias S2MktConc S2NotImplemented
typealias S2Def2 S2NotImplemented
typealias S2LifeMorb S2NotImplemented
typealias S2LifeRevision S2NotImplemented
typealias S2Health S2NotImplemented
typealias S2NonLife S2NotImplemented

"""
Projection parameters for S2-calculations without
S2 calibration parameters

  - `t_0::Int`: Year in which projection starts
  - `T::Int`: Year in which projection ends
  - `cap_mkt::CapMkt`: Capital market
  - `invs_par::Vector{Any}`:
    Parameters for construction of InvPort object with the
    exception of `t_0::Int`, `T::Int`, `cap_mkt::CapMkt`
  - `l_ins::LiabIns`:
    Liability porfolio (without other liabilites)
  - `l_other::LiabOther`: Other liabilities
  - `dyn::Dynamic`: Dynamic projection parameters
  - `tax_rate::Float64`: (constant) tax rate during projection
  - `tax_credit_0::Float64`: Initial tax credit
"""
type ProjParam
  t_0::Int
  T::Int
  cap_mkt::CapMkt
  invs_par::Vector{Any}
  l_ins::LiabIns
  l_other::LiabOther
  dyn::Dynamic
  tax_rate::Float64
  tax_credit_0::Float64
end

## solvency 2 market risk =======================================
"""
Solvency 2 market risk: interest rate risk

  - `shock_object::Symbol`: Name of object type to be shocked
  - `shock::Dict{Symbol, Any}`: Interest rate shocks
  - `spot_up_abs_min::Float64`: Minumum upwards shock
  - `balance::DataFrame`:  Shocked balance sheets
  - `scr::Vector{Float64}`: Gross and net SCR
  - `scen_up::Bool`:
    Marker whether interest was shocked upwards to obtain SCR
"""
type S2MktInt <: S2Module
  shock_object::Symbol
  shock::Dict{Symbol, Any}
  spot_up_abs_min::Float64
  balance::DataFrame
  scr::Vector{Float64}
  scen_up::Bool
end

"""
Solvency 2 market risk: equity risk

  - `shock_object::Symbol`: Name of object type to be shocked
  - `eq2type::Dict{Symbol, Symbol}`:
    Equity type for each invested stock index
  - `shock::Dict{Symbol, Any}`: Equity shocks
  - `balance::DataFrame`:  Shocked balance sheets
  - `corr::Matrix{Float64}`: Correlation matrix for equity types
  - `scr::Vector{Float64}`: Gross and net SCR
"""
type S2MktEq <: S2Module
  shock_object::Symbol
  eq2type::Dict{Symbol, Symbol}
  shock::Dict{Symbol,Any}
  balance::DataFrame
  corr::Matrix{Float64}
  scr::Vector{Float64}
end

"""
Solvency 2 market risk

  - `mds::Vector{S2Module}`: Sub-modules for market risk
  - `corr_up::Matrix{Float64}`:
    Correlation matrix if upwards shock was used in the
    calculation of the SCR for interest risk
  - `corr_down::Matrix{Float64}`:
    Correlation matrix if downwards shock was used in the
    calculation of the SCR for interest risk
  - `scr::Vector{Float64}`: Gross and net SCR
"""
type S2Mkt <: S2Module
  mds::Vector{S2Module}
  corr_up::Matrix{Float64}
  corr_down::Matrix{Float64}
  scr::Vector{Float64}
end

## solvency 2 default risk ======================================
"""
Solvency 2 default risk of type 1

  - `tlgd::Vector{Float64}`: Total loss given default per rating
  - `slgd::Vector{Float64}`:
    Squared loss given default per rating
  - `u::Matrix{Float64}`: Matrix u for variance calculation
  - `v::Vector{Float64}`: Vector v for variance calculation
  - `scr_par::Dict{Symbol, Vector{Float64}}`:
    Parameters for SCR calculation depending on whether the
    normalized standard deviation is low, medium, or high
  - `scr::Vector{Float64}`: Gross and net SCR
"""
type S2Def1 <: S2Module
  tlgd::Vector{Float64}
  slgd::Vector{Float64}
  u::Matrix{Float64}
  v::Vector{Float64}
  scr_par::Dict{Symbol, Vector{Float64}}
  scr::Vector{Float64}
end

"""
Solvency 2 default risk

  - `mds::Vector{S2Module}`: Sub-modules for default risk
  - `corr::Matrix{Float64}`: Correlation matrix for default risk
  - `scr::Vector{Float64}`: Gross and net SCR
"""
type S2Def <: S2Module
  mds::Vector{S2Module}
  corr::Matrix{Float64}
  scr::Vector{Float64}
end

## solvency 2 life risk =========================================
"""
Solvency 2 life risk: biometric risks: mortality, longevity,
surrender

  - `shock_object::Symbol`: Name of object type to be shocked
  - `shock::Dict{Symbol, Any}`: biometric shocks
  - `shock_param::Dict{Symbol, Float64}`:
    Additional parameters for calculating shocks, may be empty
  - `balance::DataFrame`:  Shocked balance sheets
  - `mp_select::Dict{Symbol, Vector{Bool}}`:
    Indicator which model points are selected to be shocked
  - `scr::Vector{Float64}`: Gross and net SCR
"""
type S2LifeBio <: S2Module
  shock_object::Symbol
  shock::Dict{Symbol, Any}
  shock_param::Dict{Symbol, Float64}
  balance::DataFrame
  mp_select::Dict{Symbol, Vector{Bool}}
  scr::Vector{Float64}
end

"""
Solvency 2 life risk: expense risk

  - `shock_object::Symbol`: Name of object type to be shocked
  - `shock::Dict{Symbol, Any}`: biometric shocks
  - `shock_param::Dict{Symbol, Float64}`:
    Additional parameters for calculating shocks
  - `balance::DataFrame`:  Shocked balance sheets
  - `scr::Vector{Float64}`: Gross and net SCR
"""
type S2LifeCost <: S2Module
  shock_object::Symbol
  shock::Dict{Symbol, Any}
  shock_param::Dict{Symbol, Float64}
  balance::DataFrame
  scr::Vector{Float64}
end

"""
Solvency 2 life risk

  - `mds::Vector{S2Module}`: Sub-modules for life risk
  - `corr::Matrix{Float64}`: Correlation matrix for life risk
  - `scr::Vector{Float64}`: Gross and net SCR
"""
type S2Life <: S2Module
  mds::Vector{S2Module}
  corr::Matrix{Float64}
  scr::Vector{Float64}
end

## solvency 2 operational risk ==================================

"""
Solvency 2 operational risk

  - `fac::Dict{Symbol, Float64}`: Factors for SCR calculation
  - `prem_earned::Float64`: Earned premium (previous year)
  - `prem_earned_prev::Float64`:
    Earned premium (year before previous year)
  - `tp::Float64`:
    technical provisions (incl. bonus but without unit linked)
  - `comp_prem::Float64`:
    Temporary result, SCR component based on technical provisions
  - `comp_tp::Float64`:
    Temporary result, SCR component based on earned premiums
  - `cost_ul::Float64`: Costs for unit linked
  - `scr::Vector{Float64}`: SCR
"""
type S2Op <: S2Module
  fac::Dict{Symbol, Float64}
  prem_earned::Float64
  prem_earned_prev::Float64
  tp::Float64
  comp_prem::Float64
  comp_tp::Float64
  cost_ul::Float64
  scr::Float64
end

## solvency 2 total =============================================
"""
Solvency 2 (total)

  - `mds::Vector{S2Module}`:
    Solvency 2 modules (with the esception of operational risk)
  - `balance::DataFrame`:  Unshocked best estimate balance sheet
  - `corr::Matrix{Float64}`:
    Correlation matrix for calculation of BSCR
  - `bscr::Vector{Float64}`: Gross and net BSCR
  - `adj_tp::Float64`:
    Adjustment for riskmitigating effect from discretionary bonus
  - `adj_dt::Float64`:
    Adjustment for riskmitigating effect from deferred taxes
  - `op::S2Op`: Solvency 2 module for operational risk
  - `scr::Float64`: SCR
  - `invest_mod::Float64`:
    Modified investments for calculation of the SCR-ratio
  - `liabs_mod::Float64`:
    Modified liabilities for calculation of the SCR-ratio
  - `scr_ratio::Float64`: SCR-ratio
"""
type S2 <: S2Module
  mds::Vector{S2Module}
  balance::DataFrame
  corr::Matrix{Float64}
  bscr::Vector{Float64}
  adj_tp::Float64
  adj_dt::Float64
  op::S2Op
  scr::Float64
  invest_mod::Float64
  liabs_mod::Float64
  scr_ratio::Float64
end
