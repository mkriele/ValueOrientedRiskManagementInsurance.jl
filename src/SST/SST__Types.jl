export Stress, Asset, Liabilities, StockIndex, ZeroBond,
       RiskFactor, SSTCapMkt

type Stress
  n::Int                ## number of scenarios
  name::Vector{String}  ## name of scenario
  target::Vector{Bool}  ## is the particular scenario included?
  prob::Vector{Float64} ## probability of the scenario
  Δx::Array{Float64,2}  ## impact of scenario
  Δrtk::Vector{Float64} ## effect on RTK
end

type SSTCapMkt             ## Capital market
  spot::Vector{Float64}
  stock_increase::Float64
end

abstract Asset

type ZeroBond <: Asset
  nom::Float64
  τ::Int
  index::Int            ## index for risk factor
end

type StockIndex <: Asset
  nom::Float64
  index::Int            ## index for risk factor
end


type Liabilities
  B_PX::Vector{Float64} ## survival benefits
  qx::Vector{Float64}
  index_spot::Vector{Int}   ## index for risk factor
  index_mort::Vector{Int}   ## index for risk factor
end

type RiskFactor
  Σ::Array{Float64,2}   ## covariance matrix for risk factors
  x0::Vector{Float64}   ## expected risk factor:  x0 = E(x)
  h::Vector{Float64}
  add::Vector{Bool}
  index_spot::Vector{Int}
  index_stock::Vector{Int}
  index_mort::Vector{Int}
end

