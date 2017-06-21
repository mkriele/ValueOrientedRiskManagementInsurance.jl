export HWBS, distr, rand_rw, zb

using Distributions
using DataFrames
using Dierckx

"""
The 1+1 Hull-White Black-Scholes model for both
real world and risk-neutral projections. The time line is the
interval [0,T], discretized into nᵀ steps of the same length.
"""
type HWBS
  "Number of steps per year"
  nʸ::Int
  "Total number of steps"
  nᵀ::Int
  "Step width"
  δt::Real
  "HW: initial interest rate"
  r₀::Real
  "HW: real world expected interest rate"
  r∞::Real
  "HW:convergence parameter"
  a::Real
  "HW: volatility parameter"
  σ::Real
  "BS: initial index price"
  S₀::Real
  "BS: expected real world yield from index"
  μ∞::Real
  "BS: volatility parameter"
  η::Real
  "BS/HW: correlation"
  ρ::Real
  "HW: spot rates etc. (risk neutral data from prices)"
  prices::DataFrame
end
"""
The input `df_spot` has at least two columns:
- `:t`: maturities for spot rates.
- `:spot`: spot rates at time zero for corresponding
    maturities in `df_spot[:t]`

The maximum time `T` is calculated as `maximum(df_spot[:t])`.
"""
function HWBS(p::Dict{Symbol, Real}, df_spot::DataFrame)
  HWBS( convert(Int, round(p[:nʸ],0) ),
        p[:r₀], p[:r∞], p[:a], p[:σ],
        p[:S₀], p[:μ∞], p[:η],
        p[:ρ],
        df_spot)
end

function HWBS(nʸ::Int,
              r₀::Real, r∞::Real, a::Real, σ::Real,
              S₀::Real, μ∞::Real, η::Real,
              ρ::Real,
              df_spot::DataFrame)
  δt = 1 / nʸ
  nᵀ = maximum(df_spot[:t]) * nʸ
  func_spot =
    Spline1D( convert(Array, df_spot[:t]),
              convert(Array, df_spot[:spot]))
  # Durations and spot rates
  prices =
    DataFrame(t = collect(1:nᵀ) * δt,
              spot = func_spot(collect(1:nᵀ) * δt))
  prices[:zb] = exp(-prices[:t] .* prices[:spot])
  # Instantaneous forward rates
  # Approximate last forward with previous forward
  prices[:forw] =
    vcat( prices[1:(end-1), :t] .*
            diff(prices[:spot]) ./ diff(prices[:t]),
          prices[end, :t] *
            diff(prices[:spot])[end] /
            diff(prices[:t])[end])
  ## the differences are not long enough and the last entry
  ## is not correct
  ## Approximate correction: repeat the last (correct) value twice
  tmp =  diff(prices[:forw]) ./ diff(prices[:t])
  prices[:dforw] =
    vcat(tmp[1:end-1], tmp[end-1], tmp[end-1])

  HWBS( nʸ, nᵀ, δt, df_spot[1,:spot], r∞, a, σ,
        S₀, μ∞, η, ρ, prices)
end

function ev_rw(cm::HWBS, t::Real)
  [(cm.r₀-cm.r∞) * exp(-cm.a * t) + cm.r∞,
  log(cm.S₀) + (cm.μ∞ - cm.η * cm.η /2)* t]
end

function cov_rw(cm::HWBS, t::Real)
  cov11 = cm.σ * cm.σ / (2cm.a) * (1 - exp(-2cm.a * t))
  cov12 = cm.ρ * cm.σ * cm.η / cm.a * (1 - exp(-cm.a * t))
  cov22 = cm.η * cm.η * t
  [cov11 cov12; cov12 cov22]
end

"""
Real world distribution of the HWBS model
- first component: ``r_t``
- second component: ``\\ln S_t``  **(ln may not be expected)**
"""
function distr_rw(cm::HWBS, t::Real)
  MvNormal(ev_rw(cm, t), cov_rw(cm, t))
end

"""
Real world random values for the HWBS model
- first component: ``r_t``
- second component: ``S_t``
"""
function rand_rw(cm::HWBS, t::Real, n)
  tmp = rand(distr_rw(cm, t), n)
  tmp[2,:] = exp.(tmp[2,:])
  tmp
end

"""
Price of zero bond which matures at time t₂, based on information at time 0.
"""
function zb(cm::HWBS, t₂::Real)
    n₂ = convert(Int, round(t₂ * cm.nʸ, 0))
  if n₂ > nrow(cm.prices)
    error("t₂ too large")
  else
    return cm.prices[n₂, :zb]
  end
end

"""
Price of zero bond which matures at time ``t₂``, as estimated  at time ``t₁``, under the condition that ``r_{t_1} = r_1``.
"""
function zb(cm::HWBS, r₁::Real, t₁::Real, t₂::Real)
  Δt = t₂-t₁
  n₁ = convert(Int, round(t₁ * cm.nʸ, 0))
  n₂ = convert(Int, round(t₂ * cm.nʸ, 0))
  if n₂ > nrow(cm.prices)
    error("t₂ too large")
  elseif Δt < 0
    error("t₂ < t₁")
  elseif n₁ == 0
    return zb(cm, t₂)
  else
    B = (1 - exp(-cm.a * Δt)) / cm.a
    A = log(cm.prices[n₂, :zb]) - log(cm.prices[n₁, :zb]) +
        cm.prices[n₁, :forw] * B -
        cm.σ^2 / (4 * cm.a) * B^2 * (1 - exp(-2 * cm.a * Δt))
    return exp(A) * exp(-B * r₁)
  end
end
