using Distributions

export value, delta, gammamatrix, Δrtk, rΔrtk, rtk, srtk, aggrstress,
UP, DOWN

## Constructors -------------------------------------------------
function RiskFactor(σ::Vector{Float64},
                    corr::Array{Float64, 2},
                    x0::Vector{Float64},
                    h::Vector{Float64},
                    add::Vector{Bool}
                    )
  return RiskFactor((σ * σ') .* corr, x0, h, add)
end

## Valuation ----------------------------------------------------
function value(t::Int,
               zb::ZeroBond,
               x::Vector{Float64},
               rf::RiskFactor,
               cap_mkt::SSTCapMkt)
  x_spot = cap_mkt.spot[zb.index] + x[zb.index]
  val = zb.nom / (1 + x_spot)^zb.τ
  if t != 0
    val *= (1 + cap_mkt.spot[t] + x[t])^t
  end
  return val
end

function value(t::Int,
               si::StockIndex,
               x::Vector{Float64},
               rf::RiskFactor,
               cap_mkt::SSTCapMkt)
  return x[si.index] * si.nom * (1 + cap_mkt.stock_increase)^t
end

function value(t::Int,
               assts::Vector{Asset},
               x::Vector{Float64},
               rf::RiskFactor,
               cap_mkt::SSTCapMkt)
  val = 0.0
  for asset in assts
    val += value(t, asset, x, rf, cap_mkt)
  end
  return val
end

function value(t::Int,
               liabs::Liabilities,
               x::Vector{Float64},
               rf::RiskFactor,
               cap_mkt::SSTCapMkt)
  T = length(cap_mkt.spot)
  x_spot = cap_mkt.spot[1:T] + x[1:T]
  x_mort = x[liabs.index_mort]
  val = 0.0
  for τ in 1:length(liabs.B_PX)
    val +=
      prod(1 .- liabs.qx[1:τ] .* x_mort) *
      liabs.B_PX[τ] / (1 + x_spot[τ])^τ
  end
  if t != 0
    val *= (1 + cap_mkt.spot[t] + x[t])^t
  end
  return val
end

rtk(t::Int,
    assets::Vector{Asset},
    liabs::Liabilities,
    x::Vector{Float64},
    rf::RiskFactor,
    cap_mkt::SSTCapMkt) =
  value(t, assets, x, rf, cap_mkt) -
  value(t,liabs, x, rf, cap_mkt)

## capital calculation ------------------------------------------
const UP, DOWN = 1, -1

"calculates (linear) sensitivities for rtk"
function srtk(shift::Int,
             assets::Vector{Asset},
             liabs::Liabilities,
             rf::RiskFactor,
             cap_mkt::SSTCapMkt)
  x = deepcopy(rf.x0)
  n = length(x)
  rtk_shift = Array(Float64, n)
  for i = 1:n
    x[i] += shift * rf.h[i]
    rtk_shift[i] = rtk(1, assets, liabs, x, rf, cap_mkt)
    x[i] -= shift * rf.h[i]  ## restore old value for y[i]
  end
  return rtk_shift
end

"calculates quadratic sensitivities for rtk"
function srtk(shift_1::Int,
             shift_2::Int,
             assets::Vector{Asset},
             liabs::Liabilities,
             rf::RiskFactor,
             cap_mkt::SSTCapMkt)
  x = deepcopy(rf.x0)
  n = length(x)
  rtk_shift_shift = Array(Float64, n, n)
  for i = 1:n
    for k = 1:n
      x[i] += shift_1 * rf.h[i]
      x[k] += shift_2 * rf.h[k]
      rtk_shift_shift[i, k] =
        rtk(1, assets, liabs, x, rf, cap_mkt)
      x[i] -= shift_1 * rf.h[i]  ## restore old value for x[i]
      x[k] -= shift_2 * rf.h[k]  ## restore old value for x[i]
    end
  end
  return rtk_shift_shift
end

Δ(rf::RiskFactor) =
  Float64[rf.add[i] ?  rf.h[i] : rf.x0[i] * rf.h[i]
          for i = 1:length(rf.x0)]

"calculates δ-vector"
delta(assets::Vector{Asset},
      liabs::Liabilities,
      rf::RiskFactor,
      cap_mkt::SSTCapMkt) =
  (srtk(UP, assets, liabs, rf, cap_mkt) -
     srtk(DOWN, assets, liabs, rf, cap_mkt)) ./ (2Δ(rf))

"calculates Γ-matrix"
function gammamatrix(assets::Vector{Asset},
               liabs::Liabilities,
               rf::RiskFactor,
               cap_mkt::SSTCapMkt)
  Δx = Δ(rf)
  rtk_uu = srtk(UP, UP, assets, liabs, rf, cap_mkt)
  rtk_ud = srtk(UP, DOWN, assets, liabs, rf, cap_mkt)
  rtk_du = srtk(DOWN, UP, assets, liabs, rf, cap_mkt)
  rtk_dd = srtk(DOWN, DOWN, assets, liabs, rf, cap_mkt)
  Γ_diag =
    (srtk(UP, assets, liabs, rf, cap_mkt) +
       srtk(DOWN, assets, liabs, rf, cap_mkt) -
       2rtk(1,assets, liabs, rf.x0, rf, cap_mkt)) ./ (Δx .* Δx)
  Γ = Array(Float64, length(rf.x0), length(rf.x0))
  for i = 1:length(rf.x0)
    for k = 1:(i-1)
      Γ[i,k] = (rtk_uu[i,k] -
                  rtk_ud[i,k] -
                  rtk_du[i,k] +
                  rtk_dd[i,k]) / (4 * Δx[i] * Δx[k])
      Γ[k,i] = Γ[i,k]
    end
    Γ[i,i] = Γ_diag[i]
  end
  return Γ
end

"calculates the shocked rtk based on shocks Δx"
Δrtk(Δx::Vector{Float64},
     δ::Vector{Float64},
     Γ::Matrix{Float64}) = (Δx ⋅ δ + 0.5 * Δx' * Γ * Δx)[1]

## random values for Δrtk
function rΔrtk(n_scen::Int,
               assets::Vector{Asset},
               liabs::Liabilities,
               rf::RiskFactor,
               cap_mkt::SSTCapMkt,
               x_index::Vector{Int},
               )
  r_Δx = rand(MvNormal(zeros(Float64, length(x_index)),
                       rf.Σ[x_index, x_index]),
              n_scen)
  r_Δrtk = Array(Float64, n_scen)
  δ = delta(assets, liabs, rf, cap_mkt)[x_index]
  Γ = gammamatrix(assets, liabs, rf, cap_mkt)[x_index, x_index]
  for mc in 1:n_scen
    r_Δrtk[mc] = Δrtk(r_Δx[:, mc], δ, Γ)
  end
  return r_Δrtk
end

##  Aggregate stress-scenarios to randomly generated values
##    r_Δrtk_no_stress.
##  We use the same approximation as in the shift method,
##    namely, two different stress scenarios cannot happen
##    in the same year.
function aggrstress(stress::Stress, r_Δrtk_no_stress)
  r_Δrtk = deepcopy(r_Δrtk_no_stress)
  n_scen = length(r_Δrtk)
  i = 0
  for scen = 1:stress.n
    if stress.target[scen]
      n_adj = floor(Integer, n_scen * stress.prob[scen])
      for j = 1:n_adj
        r_Δrtk[i + j] += min(0, stress.Δrtk[scen])
      end
      i += n_adj
    end
  end
  return r_Δrtk
end
