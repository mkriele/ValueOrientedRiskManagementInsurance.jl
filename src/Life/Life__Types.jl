export CapMkt, DetermProcess, Stock, RiskFreeRate
export Invest, InvestStock, InvestCash
export InvestGroup, IGCost, IGStock, IGCash
export InvPort
export Alloc
export Product, ModelPoint, LiabIns
export LiabOther, Debt
export Dynamic
export Projection


## capital market -----------------------------------------------

"""
A stochastic or deterministc process representing a part of the
Capital market. It always has an object `x` representing the
development in time and an initial value `x_0`.
"""
abstract Process

"""
Deterministic process, currently the only type of process
implemented

 - `x::Vector{Float64}`: process value for each Year
 - `yield_0::Float64`:
   Yield from process during year `t_0`
   (before the start of the projection)
"""
abstract DetermProcess <: Process

"""
A deterministic process representing a stock index. Additional
objects:

 - `x_0::Float64`: Initial value of process
"""
type Stock <: DetermProcess
  x_0::Float64
  x::Vector{Float64}
  yield_0::Float64
end

"A deterministic process representing a short rate"
type RiskFreeRate <: DetermProcess
  x::Vector{Float64}
  yield_0::Float64
end

"""
A collection of `Process` objects which represents a
capital market. Currently we only have two deterministic
processes, `stock::Stock` and `rfr::RiskFreeRate`.
"""
type CapMkt
  stock::Stock
  rfr::RiskFreeRate
end

## investments --------------------------------------------------

"""
Represents an individual investments and its development in
time.

  - `name::Symbol`:  Name of investment
  - `proc::Process`: Associated capital market process
  - `mv_0::Float64`: Initial market value
  - `mv`:            Development of market value in time
"""
abstract Invest             ## investment

dict_ig = Dict{Symbol, Symbol}(:IGStock => :InvestStock,
                               :IGCash => :InvestCash)

"`Invest` object associated with a `Stock` process"
type InvestStock <: Invest
  name::Symbol
  proc::Stock               ## deterministic process
  mv_0::Float64             ## initial market value
  mv::Vector{Float64}       ## market value eoy
end

"""
`Invest` object associated with a `RiskFreeRate` process.
Additional objects:

  - `lgd::Float64`: Loss given default
  - `cqs::Symbol`:  Rating of counter party
"""
type InvestCash <: Invest
  name::Symbol
  proc::RiskFreeRate        ## deterministic process
  mv_0::Float64             ## initial market value
  mv::Vector{Float64}       ## market value eoy
  lgd::Float64              ## loss given default
  cqs::Symbol               ## rating of counter-party
end

"""
A group of `Invest` objects of the same sub-type

 - `investments::Vector{Invest}`: Invest objects of the same type
 - `mv_0::Float64`: Initial total market value of invest group
 - `mv`: Development of market value of invest group in time
 - `alloc::Alloc`: Allocation of investments within invest group
 - `cost::InvestCost`: investment costs of invest group
"""
abstract InvestGroup        ## investment group

"""
Allocation of investments within an `InvestGroup` object `ig`

  - `name::Vector{Symbol}`:
    Name of investments within invest group (the indices
    correpond to the indices of `ig.investments`)
  - `total::Vector{Float64}`:
    Development of allocation to invest group in time
  - `all::Matrix{Float64}`:
    `all[τ, i]` is the allocation to `ig.investments[i]` during
    year `τ`
  - `lgd::Vector{Float64}`:
    Loss given default for each counter party (where applicable)
  - `cqs::Vector{Symbol}`:
    Rating for each counter party (where not applicable: `:na`)
"""
type Alloc  ## asset allocation for each year
  name::Vector{Symbol}      ## names of investments within group
  total::Vector{Float64}    ## alloc. for the whole invest. group
  all::Matrix{Float64}    ## alloc. for each asset within ig
  lgd::Vector{Float64}       ## loss given default for each cp
  cqs::Vector{Symbol}       ## rating for each cp (or "na")
end

"""
Investment costs for an `InvestGroup` object `ig`

  - `rel::Vector{Float64}`:
    Relative costs per year including inflation
  - `abs::Vector{Float64}`:
    Absolute costs per year including inflation
  - `infl_rel::Vector{Float64}`: Inflation for relative costs
  - `infl_abs::Vector{Float64}`: Inflation for absolute costs
  - `cum_infl_rel::Vector{Float64}`:
    Cumulative inflation for relative costs
  - `cum_infl_abs::Vector{Float64}`:
    Cumulative inflation for absolute costs
  - `total::Vector{Float64}`:
    Total costs per year (end of year)
"""
type IGCost
  rel::Vector{Float64}
  abs::Vector{Float64}
  infl_rel::Vector{Float64}
  infl_abs::Vector{Float64}
  cum_infl_rel::Vector{Float64}
  cum_infl_abs::Vector{Float64}
  total::Vector{Float64}
end

"Invest Group for `Invest` object `InvestStock`"
type IGStock <: InvestGroup
  investments::Vector{InvestStock}
  mv_0::Float64
  mv::Vector{Float64}
  alloc::Alloc
  cost::IGCost
end

"Invest Group for `Invest` object `InvestCash`"
type IGCash <: InvestGroup
  investments::Vector{InvestCash}
  mv_0::Float64
  mv::Vector{Float64}
  alloc::Alloc
  cost::IGCost
end

"""
Investment portfolio

  - `t_0::Int`:  Year in which projection starts
  - `mv_0::Float64`: Initial market value
  - `mv_boy::Vector{Float64}`:
    Market value at the Beginning of each year
  - `mv::Vector{Float64}`: Market value at the end of each year
  - `yield::Vector{Float64}`: Investment yield per year
  - `cost::Vector{Float64}`:  Investment costs per year
  - `igs::Dict{Symbol, InvestGroup}`:
    Investment groups in investment porfolio
"""
type InvPort
  t_0::Int
  mv_0::Float64
  mv_boy::Vector{Float64}
  mv::Vector{Float64}
  yield::Vector{Float64}
  cost::Vector{Float64}
  igs::Dict{Symbol, InvestGroup}
end

## insurance liabilities ----------------------------------------

"""
Life insurance product / tariff

  - `dur::Int`: Duration of product
  - `rfr::Vector{Float64}`:
    Discount rate for pricing for each year
  - `prob::DataFrame`:
    Biometric probabilities for pricing (`:qx`, `:sx`, `:px`)
  - `β::DataFrame`:
    Premium/benefit profile (`:qx`, `:sx`, `:px`, `:prem`)
  - `λ::DataFrame`: Cost profile (`:boy`, `:eoy`)
  - `prem_norm::Float64`: Normalized premium for insured sum == 1
"""
type Product
  dur::Int
  rfr::Vector{Float64}
  prob::DataFrame
  β::DataFrame
  λ::DataFrame
  prem_norm::Float64
end

"""
Model point for the liability portfolio representing a number
of insured with the same contract parameters and the same
biometric parameters

  - `n::Int`:  Number of contracts in model point
  - `t_start::Int`:  Year in which contract has been taken out
  - `dur::Int`:
    Remaining duration relative to `t_0` (projection start)
  - `prob::DataFrame`:
    Best estimate biometric probabilities (`:qx`, `:sx:`, `:px`)
  - `lx_boy::Vector{Float64}`:
    Survivors at beginning of year
  - `lx_boy_next::Float64`: The value of lx_boy for the next year
  - `β::DataFrame`:
    Conditional benefits / premium (`:qx`, `:sx`, `:px`, `:prem`)
  - `λ::DataFrame`: Conditional costs profile (`:boy`, `:eoy`)
  - `bonus_rate_hypo::Float64`:
    Hypothetical bonus rate communicated when contract was sold
  - `rfr_price_0::Float64`: Discount rate for pricing at `t_0`
  - `rfr_price::Vector{Float64}`: Discount rate for pricing
  - `tpg_price_0::Float64`:
    Initial technical provisions for pricing
  - `tpg_price::Vector{Float64}`: Technical provisions for pricing
  - `gc::Vector{Float64}`:
    Going concern factor (fraction of policy holders per year)
  - `pension_contract::Bool`:
    Marker whether it is a pension contract (for S2-calculation)

**Time model:**

```
  project time τ:          0
  real time t:    t_start  t_0
  product time s: 0        s_0
                  |--------|-------------------|--------------
                           \-------------------/
                                     dur
                  \----------------------------/
                           product.dur
```
"""
type ModelPoint
  n::Int
  t_start::Int
  dur::Int
  prob::DataFrame
  lx_boy::Vector{Float64}
  lx_boy_next::Float64
  β::DataFrame
  λ::DataFrame
  bonus_rate_hypo::Float64
  rfr_price_0::Float64
  rfr_price::Vector{Float64}
  tpg_price_0::Float64
  tpg_price::Vector{Float64}
  gc::Vector{Float64}
  pension_contract::Bool
end

"""
Liability portfolio

  - `n::Int`: Number of model points
  - `t_0::Int`: Year in which projection starts
  - `dur::Int`: Maximum remaining duration in liability portfolio
  - `mps::Vector{ModelPoint}`: Model points in portfolio
  - `gc::Vector{Float64}`: Scaling for going concern modeling
  - `Δgc::Vector{Float64}`: `Δgc[t] = gc[t + 1] - gc[t]`
"""
type LiabIns
  n::Int
  t_0::Int
  dur::Int
  mps::Array{ModelPoint, 1}
  gc::Vector{Float64}
  Δgc::Vector{Float64}
end

## other liabilities --------------------------------------------

"""
Debt

  - `name::Symbol`: Name of this loan
  - `t_init::Int`:  Time at which loan has been taken
  - `t_mat::Int`:   Time at which loan matures
  - `τ_init::Int`:  Time at loan has been taken relative to `t_0`
  - `τ_mat::Int`:   Remaining duration of loan relative to `t_0`
  - `nominal::Float64`: Nominal amount
  - `coupon::Float64`:  Yearly coupon payment (absolute value)
"""
type Debt
  name::Symbol
  t_init::Int
  t_mat::Int
  τ_init::Int
  τ_mat::Int
  nominal::Float64
  coupon::Float64
end

"""
  Portfolio of other liabilities

  - `subord::Vector{Debt}`:  Vector with subordinated debt
"""
type LiabOther
  subord::Vector{Debt}
end

## dynamics -----------------------------------------------------

"""
Dynamic projection parameters for a `Projection` object

  - `bonus_factor::Float64`:  Bonus factor
  - `quota_surp::Float64`:    Desired surplus ratio
  - `free_surp_boy::Vector{Float64}`:
    Free surplus for bonus calculation
"""
type Dynamic
  bonus_factor::Float64
  quota_surp::Float64
  free_surp_boy::Vector{Float64}
end

## cashflow projection ------------------------------------------
"""
  Projection of a life insurer

  - `t_0::Int`: Start of projection
  - `dur::Int`: Length of projection
  - `cf::DataFrame`:
    Cashflows (`:qx`,`:sx`,`:px`,`:prem`,`:λ_boy`,`:λ_eoy`,
    `:bonus`, `:invest`,`:new_debt`,`:l_other`,`:profit`,
    `:tax`, `:divid`, `:gc`, `:Δtpg`, `:cost_prov`)
  - `val_0::DataFrame`:
    Best estimate initial valuation
    (`:invest`, `:tpg`, `:l_other`, `:surplus`, `:bonus`,
    `:cost_prov`)
  - `val::DataFrame`:
    Best estimate valuation for each year
    (`:invest`, `:tpg`, `:l_other`, `:surplus`, `:bonus`,
    `:cost_prov`)
  - `tax_rate::Float64`:  Rate at which profit is taxed
  - `tax_credit_0::Float64`: Initial tax credit
  - `tax_credit::Vector{Float64}`: Tax credit for each year
  - `fixed_cost_gc::Vector{Float64}`:
    Weighted fixed costs for going concern modeling
"""
type Projection
  t_0::Int
  dur::Int
  cf::DataFrame
  val_0::DataFrame
  val::DataFrame
  tax_rate::Float64
  tax_credit_0::Float64
  tax_credit::Vector{Float64}
  fixed_cost_gc::Vector{Float64}
end
