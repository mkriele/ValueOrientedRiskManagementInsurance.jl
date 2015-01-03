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
abstract DetermProcess
type Stock <: DetermProcess
  x_0::Float64               ## initial process value
  x::Vector{Float64}         ## process values eoy
  yield_0::Float64             ## yield of year t_0
end

type RiskFreeRate <: DetermProcess
  x::Vector{Float64}         ## process values boy
  yield_0::Float64             ## yield of year t_0
end

## our capital market only consists of one stock and the
## risk free interest rate.
type CapMkt
  stock::Stock              ## stock in the market
  rfr::RiskFreeRate         ## "short" means duration = one year
end

## investments --------------------------------------------------
## individual investments
abstract Invest             ## investment

dict_ig = Dict{Symbol, Symbol}([:IGStock => :InvestStock,
                                :IGCash => :InvestCash])

type InvestStock <: Invest
  name::Symbol
  proc::Stock               ## deterministic process
  mv_0::Float64             ## initial market value
  mv::Vector{Float64}       ## market value eoy
end

type InvestCash <: Invest
  name::Symbol
  proc::RiskFreeRate        ## deterministic process
  mv_0::Float64             ## initial market value
  mv::Vector{Float64}       ## market value eoy
  lgd::Float64              ## loss given default
  cqs::Symbol               ## rating of counter-party
end

## investments per investment group
abstract InvestGroup        ## investment group

type Alloc  ## asset allocation for each year
  name::Vector{Symbol}      ## names of investments within group
  total::Vector{Float64}    ## alloc. for the whole invest. group
  all::Array{Float64, 2}    ## alloc. for each asset within ig
  lgd::Vector{Float64}       ## loss given default for each cp
  cqs::Vector{Symbol}       ## rating for each cp (or "na")
end

type IGCost
  rel::Vector{Float64}      ## rel. costs per year incl. infl.
  abs::Vector{Float64}      ## abs. costs per period incl. infl.
  total::Vector{Float64}    ## total costs per year eoy
end

type IGStock <: InvestGroup ## investments in stock indices
  investments::Vector{InvestStock}
  mv_0::Float64             ## inital mv eoy: stock investm.
  mv::Vector{Float64}       ## mv eoy: stock investm.
  alloc::Alloc              ## investment allocation
  cost::IGCost              ## investment costs
end

type IGCash <: InvestGroup  ## investments in stock indices
  investments::Vector{InvestCash}
  mv_0::Float64             ## inital mv eoy: cash investm.
  mv::Vector{Float64}       ## mv eoy: cash investm.
  alloc::Alloc              ## investment allocation
  cost::IGCost              ## investment costs
end

## All investments
type InvPort               ## all investments
  t_0::Int                  ## time at which projection starts
  mv_0::Float64             ## initial market value
  mv_boy::Vector{Float64}   ## market value boy
  mv::Vector{Float64}       ## market value eoy
  yield::Vector{Float64}    ## investment yield per year
  cost::Vector{Float64}     ## investment costs
  igs::Dict{Symbol, InvestGroup}  ## investment groups
end

## insurance liabilities ----------------------------------------
type Product
  dur::Int                  ## duration of product
  rfr::Vector{Float64}      ## discount rate
  prob::DataFrame           ## biometric probabilities (pricing)
  β::DataFrame              ## benefit profile
  λ::DataFrame              ## cost profile
  prem_norm::Float64        ## norm. premium for insured sum = 1
  #   tp_norm::Vector{Float}    ## norm. techn. provisions (pricing)
end

type ModelPoint             ## liability model point
  n::Int                    ## number of contracts in model point
  dur::Int                  ## remaining duration
  prob::DataFrame           ## biometric probabilities (be)
  lx_boy::Vector{Float64}   ## survivors at boy
  lx_boy_next::Float64      ## lx_boy for next year
  β::DataFrame              ## conditional benefits / Premium
  λ::DataFrame              ## conditional costs profile
  bonus_rate_hypo::Float64  ## hypothetical bonus rate
  rfr_price::Vector{Float64} ## discount rate for pricing
  tpg_price_0::Float64       ## initial techn. prov. for pricing
  tpg_price::Vector{Float64} ## technical provisions for pricing
  gc::Vector{Float64}       ## going concern (fraction per year)
end

type LiabIns               ## liability portfolio
  n::Int                    ## number of model points
  t_0::Int                  ## time at which portfolio is defined
  dur::Int                  ## maximum duration in Liab_port
  mps::Array{ModelPoint, 1} ## model points in portfolio
  gc::Vector{Float64}       ## going concern scaling
end

## other liabilities --------------------------------------------
type Debt
  name::Symbol
  t_0::Int          ## time at which loan has been taken
  t_mat::Int        ## time at which loan matures
  τ_0::Int          ## time at loan has been taken rel t_0
  τ_mat::Int               ## duration of loan relative to t_0
  nominal::Float64       ## nominal amount
  coupon::Float64        ## coupon payment (absolute value)
end

type LiabOther           ## other liabilities
  subord::Vector{Debt}   ## subordinated debt
end

## dynamics -----------------------------------------------------
type Dynamic
  bonus_factor::Float64
  divid_factor::Float64
  surp_factor::Float64
  surp_0::Float64            ## initial surplus for bonus calc.
  surp_threshold_0::Float64  ## init. surplus threshold for bonus
  surp::Vector{Float64}      ## surplus threshold for bonus calc.
  surp_threshold::Vector{Float64}  ## surplus threshold for bonus
end

## cashflow projection ------------------------------------------
type Projection
  t_0::Int                 ## start of projection
  dur::Int                 ## length of projection
  cf::DataFrame            ## cashflows
  val_0::DataFrame         ## initial valuation (best estimate)
  val::DataFrame           ## valuation (best estimate)
  tax_rate::Float64        ## tax rate for profit
  tax_credit_0::Float64
  tax_credit::Vector{Float64}
  fixed_cost_gc::Vector{Float64}  ## weighted fixed costs
end

