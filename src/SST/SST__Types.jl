export Stress, Asset, Liabilities, StockIndex, ZeroBond,
       RiskFactor, SSTCapMkt

"""
Stress scenario

  - `n::Int`: Number of scenarios
  - `name::Vector{AbstractString}`: Name of scenario
  - `target::Vector{Bool}`: Is the particular scenario included?
  - `prob::Vector{Float64}`: Probability of the scenario
  - `Δx::Array{Float64,2}`: Impact of scenario
  - `Δrtk::Vector{Float64}`: Effect on RTK
"""
type Stress
  n::Int
  name::Vector{AbstractString}
  target::Vector{Bool}
  prob::Vector{Float64}
  Δx::Array{Float64,2}
  Δrtk::Vector{Float64}
end

"""
Capital market

  - `spot::Vector{Float64}`: Risk-free spot rate curve
  - `stock_increase::Float64`: Expected relative increase
     of stocks index per year
"""
type SSTCapMkt
  spot::Vector{Float64}
  stock_increase::Float64
end

abstract Asset

"""
Zero bond (an instance of an asset)

  - `nom::Float64`: Nominal value of the zero bond
  - `τ::Int`: Time to maturity (full years)
  - `index::Int`: Index of this risk factor
"""
type ZeroBond <: Asset
  nom::Float64
  τ::Int
  index::Int
end

"""
Investment in stock index

  - `nom::Float64`: Initial value of investment
  - `index::Int`: Index of this risk factor
"""
type StockIndex <: Asset
  nom::Float64
  index::Int
end

"""
Liability portfolio consisting of a single insurance contract

  - `B_PX::Vector{Float64}`: Survival benefits
  - `qx::Vector{Float64}`: Mortality probabilities
  - `index_spot::Vector{Int}`:
     Indices of spot rates used for discounting
  - `index_mort::Vector{Int}`: Index of this risk factor
"""
type Liabilities
  B_PX::Vector{Float64}
  qx::Vector{Float64}
#   index_spot::Vector{Int}
  index_mort::Vector{Int}
end

"""
Risk factor()

  - `Σ::Array{Float64,2}`: Covariance matrix for risk factors
  - `x0::Vector{Float64}`: Expected value of risk factor:
     x0 = E(x)
  - `h::Vector{Float64}`: Sensitivity for each risk factor to
     calculate derivatives of RTK
  - `add::Vector{Bool}`:
     Is the sensitivity additive or multiplicative?
  - `index_spot::Vector{Int}`:
  - `index_stock::Vector{Int}`:
  - `index_mort::Vector{Int}`:
"""
type RiskFactor
  Σ::Array{Float64,2}
  x0::Vector{Float64}
  h::Vector{Float64}
  add::Vector{Bool}
#   index_spot::Vector{Int}
#   index_stock::Vector{Int}
#   index_mort::Vector{Int}
end

