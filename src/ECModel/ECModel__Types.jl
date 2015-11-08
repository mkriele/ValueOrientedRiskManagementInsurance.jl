export GaussCopula, ProfitLoss, PLInsurance, PLInvestments,
PLTotal, BusinessUnit, BuInsurance, BuInvestments, Total

type GaussCopula
  n::Int               ## number marginal distributions
  Σ::Array{Float64,2}  ## correlation matrix
  marginals::Vector{ContinuousUnivariateDistribution}
                       ## marginal distributions
end

abstract ProfitLoss

type PLInsurance <: ProfitLoss
  premium::Real
  costs::Real
  ceded::Real
  profit::Vector{Real}
  profit_mean::Real
  eco_cap::Real
  rorac::Real
end

type PLInvestments <: ProfitLoss
  invest_bop::Real
  costs::Real
  profit::Vector{Real}
  profit_mean::Real
  eco_cap::Real
  rorac::Real
end

type PLTotal <: ProfitLoss
  profit::Vector{Real}
  profit_mean::Real
  eco_cap::Real
  rorac::Real
end

abstract BusinessUnit

type BuInsurance <: BusinessUnit
  name::AbstractString
  gross::PLInsurance
  net::PLInsurance
end

type BuInvestments <: BusinessUnit
  name::AbstractString
  init::Real
  gross::PLInvestments
  net::PLInvestments
end

type Total
  gross::PLTotal
  net::PLTotal
end
