using ValueOrientedRiskManagementInsurance
using DataFrames
using Base.Test
using Distributions

println("start NonLife.jl...")

function logpar2statpar(m_s::Vector)
  Real[exp(m_s[1] + m_s[2]^2 / 2),
       √(exp(m_s[2]^2)-1) * exp(m_s[1] + m_s[2]^2 / 2)]
end

function statpar2logpar(μ_σ::Vector)
  vc = μ_σ[2]/μ_σ[1]
  Real[log(μ_σ[1]) - 0.5 * log(1+vc^2), √(log(1+vc^2))]
end

"""
   dr = ( θ(t)-ar ) dt + σ dW
"""

type HullWhite
  Δt::Real
  T::Real
  r₀::Real
  θ::Vector{Real}
  a::Real
  σ::Real
end

function HullWhite(Δt::Real,
                   r₀::Real,
                   θ::Vector{Real},
                   a::Real,
                   σ::Real)
  T = length(θ) * Δt
  HullWhite(Δt, T, r₀, θ, a, σ)
end

periods_int(t::Real, hw::HullWhite) = int(div(t, hw.Δt))
periods_frac(t::Real, hw::HullWhite) = t/hw.Δt - period(t, hw)

function r(t::Real, hw::HullWhite)
  k = periods_int(t, hw)
  μ =
    exp(-a*t)/a *
    (a*r₀ - θ[1] +
     sum([exp(-i*a*hw.Δt) * (θ[i-1]-θ[i]) for i in 2:(k+1)])) +
     exp(-a*t) * θ[k]
  std = σ * √( (1-2exp(-2a*t)) / (2a) )
  Normal(μ, std)
end

nllobs = [:fire, :liab, :theft]
n_nl = length(nllobs)
nl_names = [ucfirst(string(nllobs[i])) for i in 1:n_nl]
df_claims = Vector(n_nl)
claims = Vector(n_nl)
res = Vector(n_nl)
β = Vector(n_nl)

claimpath = "test/NonLife_Input_Claims_"
for i = 1:n_nl
  df_claims[i] =
    readtable(claimpath * nl_names[i]  * ".csv",
              header = false)
  res[i] = Mack(df_claims[i])
  β[i] = res[i].futureclaims / sum(res[i].futureclaims)
end


#
# mean, sde -> logmean, logsigma
# R[i] = LogNormal(log_mean[i], log_sigma[i])
# v⋅mack[i].β * R[i]
#
# v: stochastic discont
# β
# R[i] = stochastic total undiscounted future paiments,
#

println("...end NonLife.jl")
