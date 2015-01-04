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

# typealias S2Life S2NotImplemented
typealias S2Health S2NotImplemented
typealias S2NonLife S2NotImplemented

typealias S2MktProp S2NotImplemented
typealias S2MktSpread S2NotImplemented
typealias S2MktFx S2NotImplemented
typealias S2MktConc S2NotImplemented

typealias S2Def2 S2NotImplemented

typealias S2LifePx S2NotImplemented
typealias S2LifeMorb S2NotImplemented
typealias S2LifeSx S2NotImplemented
typealias S2LifeCost S2NotImplemented
typealias S2LifeRevision S2NotImplemented
typealias S2LifeCat S2NotImplemented


## projection parameters for s2-calculations apart from
## s2 calibration paramaters
type ProjParam
  t_0::Int
  T::Int
  cap_mkt::CapMkt
  invs_par::Array{Any, 1}
  l_ins::LiabIns
  l_other::LiabOther
  dyn::Dynamic
  proj_par::Array{Any, 1}
end

## solvency 2 market risk =======================================

type S2MktInt <: S2Module
  shock_object::Symbol            ## object type to be shocked
  shock_type::Vector{Symbol}      ## type of shock
  shock::Dict{Symbol, Any}        ## spot_up, spot_down,
  spot_up_abs_min::Float64        ## minimum upwards shock
  balance::DataFrame              ## shocked balance sheets
  #   corr::Matrix{Float64}
  scr::Vector{Float64}
  scen_up::Bool                   ## interest shock up => SCR
end

type S2MktEq <: S2Module
  shock_object::Symbol            ## object type to be shocked
  shock_type::Vector{Symbol}      ## type of shock
  eq2type::Dict{Symbol, Symbol}   ## stock invest. => shock_type
  shock::Dict{Symbol,Any}         ## shocks[i] is always Float64,
  balance::DataFrame              ## shocked balance sheets
  corr::Matrix{Float64}
  scr::Vector{Float64}
end

type S2Mkt <: S2Module
  mds::Vector{S2Module}
  corr_up::Matrix{Float64}
  corr_down::Matrix{Float64}
  scr::Vector{Float64}
end
## solvency 2 default risk ======================================
type S2Def1 <: S2Module
  tlgd::Vector{Float64}
  slgd::Vector{Float64}
  u::Matrix{Float64}
  v::Vector{Float64}
  scr_par::Dict{Symbol, Vector{Float64}}
  scr::Vector{Float64}
end

type S2Def <: S2Module
  mds::Vector{S2Module}
  corr::Matrix{Float64}
  scr::Vector{Float64}
end

## solvency 2 life risk =========================================

type S2LifeQx <: S2Module
  shock_object::Symbol       ## object type to be shocked
  shock_type::Vector{Symbol} ## type of shock: ([:qx])
  shock::Float64
  balance::DataFrame         ## shocked balance sheets
  corr::Matrix{Float64}      ## trivial:  1
  mp_select::Vector{Bool}    ## model points selected for shock
  scr::Vector{Float64}
end

type S2Life <: S2Module
  mds::Vector{S2Module}
  corr::Matrix{Float64}
  scr::Vector{Float64}
end

## solvency 2 operational risk ==================================

type S2Op <: S2Module
  prem_earned::Float64
  prem_earned_prev::Float64
  tp::Float64
  cost_ul::Float64
  scr::Float64
end

## solvency 2 total =============================================
type S2 <: S2Module
  mds::Vector{S2Module}
  balance::DataFrame
  corr::Matrix{Float64}
  bscr::Vector{Float64}
  adj_tp::Float64
  adj_dt::Float64
  op::S2Op
  scr::Float64
end
