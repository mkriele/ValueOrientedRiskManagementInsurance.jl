export forw2spot, spot2forw

## capital market -----------------------------------------------
"""
`yield(τ, stock::Stock)`

calculates the yield of `stock` during year `τ`
"""
function yield(τ, stock::Stock)
  if τ > 1
    return stock.x[τ] / stock.x[τ - 1] - 1
  elseif τ == 1
    return stock.x[τ] / stock.x_0 - 1
  else
    return stock.yield_0
  end
end

"""
`yield(τ, rfr::RiskFreeRate)`

calculates the yield of `rfr` during year `τ`
"""
function yield(τ, rfr::RiskFreeRate)
  if τ >= 1
    return rfr.x[τ]
  else
    return rfr.yield_0
  end
end

"""
`forw2spot(f::Vector{Float64})`

calculates the spot rate `s` rate from the forward rate `f`:

  `(1+f[1])(1+f[2])...(1+f[n]) = (1+s[n])^n`
"""
forw2spot(f::Vector{Float64}) =
  cumprod(1 .+ f) .^ (1 ./ collect(1:length(f))) - 1

"""
`spot2forw(s::Vector{Float64})`

calculates the forward rate `f` rate from the spot rate `s`:

  `(1 + s[n-1])^(n-1) * (1+f[n]) = (1+s[n])^n`
"""
function spot2forw(s::Vector{Float64})
  f = zeros(Float64, length(s))
  for n in length(s):-1:2
    f[n] = (1 + s[n])^n / (1 + s[n-1])^(n - 1) -1
  end
  f[1] = s[1]
  return f
end

## investments --------------------------------------------------
"""
`project!(τ::Int, mv_boy::Float64, invest::Invest)`

project `invest::Invest` one year for a given
initial market value.

**Changed**:  `invest`
"""
function project!(τ::Int, mv_boy::Float64, invest::Invest)
  invest.mv[τ] = (1 + yield(τ, invest.proc)) * mv_boy
end

"""
`project!(τ::Int, mv_bop_total::Float64, ig::InvestGroup)`

project `ig::InvestGroup` one year for a given
initial market value.

**Changed**:  `ig`
"""
function project!(τ::Int, mv_bop_total::Float64, ig::InvestGroup)
  mv_bop = ig.alloc.total[τ] * mv_bop_total
  ig.mv[τ] = 0
  for (i, ig_invest) in enumerate(ig.investments)
    project!(τ, ig.alloc.all[τ, i] * mv_bop, ig_invest)
    ig.mv[τ] += ig_invest.mv[τ]
  end
  ig.cost.total[τ] =
    ig.cost.abs[τ] * ig.cost.cum_infl_abs[τ] +
    mv_bop * ig.cost.rel[τ] * ig.cost.cum_infl_rel[τ]
end

"""
`alloc!(τ, cap_mkt::CapMkt, invs::InvPort)`

Dynamically re-allocate investments `invs::InvPort` at the
beginning of year `τ`

**Changed**:  `invs`
"""
function alloc!(τ, cap_mkt::CapMkt, invs::InvPort)
  if τ > 1
    invs.igs[:IGStock].alloc.total[τ] =
      0.5 * (1 - exp( -max(0, dynstateavg(τ, cap_mkt))))
    invs.igs[:IGCash].alloc.total[τ] =
      1 - invs.igs[:IGStock].alloc.total[τ]
    ## we leave the allocations within each group unchanged:
    for symb in [:IGCash, :IGStock]
      for i = 1:size(invs.igs[symb].alloc.all, 2)
        invs.igs[symb].alloc.all[τ, i] =
          invs.igs[symb].alloc.all[τ-1, i]
      end
    end
  end
end

"""
`project!(τ::Int, mv_boy::Float64, invs::InvPort)`

Project `invs::InvPort`  one year for a given
initial market value.

**Changed**:  `invs`
"""
function project!(τ::Int, mv_boy::Float64, invs::InvPort)
  invs.mv[τ] = 0.0
  invs.cost[τ] = 0.0
  invs.mv_boy[τ] = mv_boy
  for ig in values(invs.igs)
    project!(τ, mv_boy, ig)
    invs.mv[τ] += ig.mv[τ]
    invs.cost[τ] += ig.cost.total[τ]
  end
  invs.yield[τ] = invs.mv[τ] / mv_boy - 1
end

## insurance liabilities  ---------------------------------------
"""
`premium(ins_sum, rfr, prob, β, λ)`

Calculate the premium of a product
"""
function premium(ins_sum, rfr, prob, β, λ)
  lx_boy = [1; cumprod(prob[:px])[1:end-1]]
  v_eoy = 1 ./ cumprod(1 .+ rfr)
  v_boy = [1; v_eoy[1:end-1]]
  num =
    sum(lx_boy .* ins_sum .*
        (v_boy .* λ[:boy] .* λ[:cum_infl] ./ (1 + λ[:infl]) +
           v_eoy .* (λ[:eoy] .* λ[:cum_infl] +
                       prob[:px] .* β[:px] +
                       prob[:qx] .* β[:qx])
         ))
  denom =
    sum(lx_boy .* β[:prem] .*
        (v_boy - v_eoy .* prob[:sx] .* β[:sx]))
  return num / denom
end

"""
`tpgrec(τ, tpg, rfr, prob, β, λ)`

One year backward recursion formula for technical provisions
of guaranteed benefits
"""
function tpgrec(τ, tpg, rfr, prob, β, λ)
  disc = 1/(1 + rfr[τ + 1])
  -β[τ + 1, :prem] +
    λ[τ + 1, :boy] *
    λ[τ + 1, :cum_infl] / (1 + λ[τ + 1, :infl]) +
    disc * (λ[τ + 1, :eoy] * λ[τ + 1, :cum_infl] +
              prob[τ + 1, :qx] * β[τ + 1, :qx] +
              prob[τ + 1, :sx] * β[τ + 1, :sx] +
              prob[τ + 1, :px] * (β[τ + 1, :px] + tpg))
end

"""
`tpg(τ, rfr, prob, β, λ)`

Technical provisions of guaranteed benefits
at the end of year `τ`
"""
function tpg(τ, rfr, prob, β, λ)
  dur = nrow(β)
  res = 0.0
  if τ >= dur
    return 0.0
  else
    for s in (dur-1):-1:τ
      res = tpgrec(s, res, rfr, prob, β, λ)
    end
    return res
  end
end

"""
`tpg(τ, rfr, mp)`

Best estimate guaranteed technical provisions
of a model point `mp` at the end of year `τ`
"""
tpg(τ, rfr, mp) = tpg(τ, rfr, mp.prob, mp.β, mp.λ)

"""
`tpgfixed(τ, rfr, fixed_cost_gc::Vector)`

Technical provisions at the end of year `τ`
for future fixed (going concern) costs
"""
function tpgfixed(τ, rfr, fixed_cost_gc::Vector)
  dur = length(fixed_cost_gc)
  if dur < τ + 1
    return 0.0
  else
    disc = cumprod(1 ./ (1 .+ rfr[(τ + 1):dur]))
    return fixed_cost_gc[(τ + 1):dur] ⋅ disc
  end
end

## other liabilities --------------------------------------------
"""
`pv(τ::Int, cap_mkt::CapMkt, debt::Debt)`

Present value of `debt` at the end of year `τ`.
Any servicing of this debt during year `τ` has occured beforehand
"""
function pv(τ::Int, cap_mkt::CapMkt, debt::Debt)
  ## calculate pv at the end of the year after servicing debt
  if (debt.τ_init > τ) | (debt.τ_mat <= τ)
    return 0.0
  else
    p_v = debt.nominal
    for s = (debt.τ_mat - 1) : -1 : τ
      p_v = (debt.coupon + p_v) / (1 + cap_mkt.rfr.x[s + 1])
    end
    return p_v
  end
end

"""
`pv(τ::Int, cap_mkt::CapMkt, l_other::LiabOther)`

present value of other liabilities `l_other` at the end of year
`τ`. Any servicing of debt within this portfolio during year `τ`
has occured beforehand
"""
function pv(τ::Int, cap_mkt::CapMkt, l_other::LiabOther)
  p_v = 0.0
  for debt in l_other.subord
    p_v += pv(τ, cap_mkt, debt)
  end
  return p_v
end

"""
`paycoupon(τ::Int, debt::Debt)`

Coupon payment for `debt` at the end of year `τ`
"""
paycoupon(τ::Int, debt::Debt) =
  (debt.τ_init <= τ <= debt.τ_mat ? debt.coupon : 0.0)

"""
`paycoupon(τ::Int, l_other::LiabOther)`

Coupon payment for all debts within the portfolio of
other liabilities `l_other` at the end of year `τ`
"""
function paycoupon(τ::Int, l_other::LiabOther)
  pay = 0.0
  for debt in l_other.subord
    pay += paycoupon(τ, debt)
  end
  return pay
end

"""
`payprincipal(τ::Int, debt::Debt)`

Payment of the principal of `debt` at the end of year `τ`,
if the debt matures at this point in time
"""
payprincipal(τ::Int, debt::Debt) =
  (τ == debt.τ_mat ? debt.nominal : 0.0)

"""
`payprincipal(τ::Int, l_other::LiabOther)`

Payment of the total principal of all debts within the
portfolio of other liabilities `l_other`, which mature at
the end of year `τ`
"""
function payprincipal(τ::Int, l_other::LiabOther)
  pay = 0.0
  for debt in l_other.subord
    pay += payprincipal(τ, debt)
  end
  return pay
end

"""
`getloan(τ::Int, debt::Debt)`

Get the nominal of `debt` at time `τ`, if it has been
taken out at this point in time, otherwise return 0.0
"""
getloan(τ::Int, debt::Debt) =
  (τ == debt.τ_init ? debt.nominal : 0.0)

"""
`getloan(τ::Int, l_other::LiabOther)`

Get the total nominal of all debts within the
portfolio of other liabilities `l_other`,
which are taken out at the beginning of year `τ`
"""
function getloan(τ::Int, l_other::LiabOther)
  nominal = 0.0
  for debt in l_other.subord
    nominal += getloan(τ, debt)
  end
  return nominal
end

"""
`goingconcern(debts::Vector{Debt}, Δgc::Vector{Float64})`

Calculates a vector of debts from an existing vector of debts
according to the going concern assumption. The total initial debt
is the same, but the new debt vectors mature earlier so that the
total nominal decreases according to the going concern factors.
We do not input the factors `gc` directly but their year on year
differences `Δgc`.
"""
function goingconcern(debts::Vector{Debt}, Δgc::Vector{Float64})
  new_debt_vec = Array(Debt, 0)
  for debt in debts
    if debt.nominal > 0.0
      τ_init = max(1, debt.τ_init)
      diff_nom = -Δgc * debt.nominal
      for τ = τ_init:debt.τ_mat
        t = debt.t_init + τ - debt.τ_init
        push!(new_debt_vec,
              Debt(debt.name,
                   debt.t_init,
                   t,
                   debt.τ_init,
                   τ,
                   diff_nom[τ],
                   debt.coupon * diff_nom[τ] / debt.nominal))
      end
    end
  end
  return(new_debt_vec)
end

"""
`goingconcern!(l_other::LiabOther, Δgc::Vector{Float64})`

Transforms a portfolio of other liabilities `l_other::LiabOther`
according to the going concern assumption. We do not input the
factors `gc` directly but their year on year differences `Δgc`.

**Changed**: `l_other`
"""
function goingconcern!(l_other::LiabOther, Δgc::Vector{Float64})
  l_other.subord = goingconcern(l_other.subord, Δgc)
end

## dynamics -----------------------------------------------------
"""
`dynstate(τ, cap_mkt::CapMkt)`

indicator for the state of the economy at the end of year `τ`
"""
dynstate(τ, cap_mkt::CapMkt) =
  yield(τ, cap_mkt.stock) / max(yield(τ, cap_mkt.rfr), eps()) - 1

"""
`dynstateavg(τ, cap_mkt::CapMkt)`

Two year average of the indicator for the state of the economy
at the end of year `τ`.
"""
dynstateavg(τ, cap_mkt::CapMkt) =
  0.5 * (yield(τ - 1, cap_mkt.stock) /
           max(yield(τ-1, cap_mkt.rfr), eps())
         + yield(τ, cap_mkt.stock) /
           max(yield(τ, cap_mkt.rfr), eps())) - 1


"""
`bonusrate(yield_eoy, rfr_price, bonus_factor)`

Helper function for dynamic bonus rate declaration
"""
bonusrate(yield_eoy, rfr_price, bonus_factor) =
  max(bonus_factor * (yield_eoy - rfr_price), 0.0)

"""
`bonusrate(τ, yield_eoy, mp::ModelPoint, dyn)`

Dynamic bonus rate declaration for a model point at the end
of year `τ`
"""
bonusrate(τ, yield_eoy, mp::ModelPoint, dyn) =
  bonusrate(yield_eoy,
            ( τ == 0 ? mp.rfr_price_0 : mp.rfr_price[τ]),
            dyn.bonus_factor)

"""
`biquotient(τ, yield_eoy, cap_mkt, invs, mp, dyn)`

Indicator for bonus rate expectation at the end of year `τ`
"""
function biquotient(τ, yield_eoy, cap_mkt, invs, mp, dyn)
  if τ ≤ mp.dur
    ind_bonus =
      yield(τ, cap_mkt.stock) /
      max(0.0,
          bonusrate(τ - 1, yield_eoy, mp, dyn) + mp.rfr_price[τ])
    ind_bonus_hypo =
      yield(0, cap_mkt.stock) /
      max(eps(), mp.bonus_rate_hypo + mp.rfr_price_0)
    return ind_bonus / ind_bonus_hypo
  else
    return 0.0
  end
end

"""
`δsx(τ, cap_mkt, invs, mp, dyn)`

Dynamic lapse probability factor to adjust the initial estimate
"""
function δsx(τ, cap_mkt, invs, mp, dyn)
  yield_eoy =
    invs.igs[:IGStock].alloc.total[τ] * yield(τ, cap_mkt.stock) +
    invs.igs[:IGCash].alloc.total[τ] * yield(τ, cap_mkt.rfr)
  if τ - 1 > mp.t_start
    bi_quot =  biquotient(τ, yield_eoy, cap_mkt, invs, mp, dyn)
  else
    bi_quot = 1
  end
  state_quot = dynstate(τ, cap_mkt) / dynstate(0, cap_mkt)
  δ_SX = 1.0
  if state_quot < 0.5
    δ_SX += 0.15
  elseif  state_quot > 2.0
    δ_SX -= 0.15
  end
  δ_SX += 0.25 * min(4.0, max(0.0, bi_quot - 1.2))
  return δ_SX
end

"""
`freesurp(dyn, invest_pre, liab)`

Helper function for the free surplus calculation
"""
freesurp(dyn, invest_pre, liab) =
 max(0, invest_pre - (1 + dyn.quota_surp) * liab)

"""
`freesurp(τ, proj::Projection, dyn)`

Free surplus for the dynamic dividend declaration
"""
freesurp(τ, proj::Projection, dyn) =
  if τ == 1
      freesurp(dyn,
               proj.val_0[1, :invest],
               proj.val_0[1, :tpg] + proj.val_0[1, :l_other])
  else
      freesurp(dyn,
               proj.val[τ-1, :invest],
               proj.val[τ-1, :tpg] + proj.val[τ-1, :l_other])
  end

"""
`update!(τ, proj::Projection, dyn::Dynamic)`

Update dynamic parameters
"""
function update!(τ, proj::Projection, dyn::Dynamic)
  if τ == 1
    dyn.free_surp_boy[τ] =
      freesurp(dyn,
               proj.val_0[1, :invest],
               proj.val_0[1, :tpg] + proj.val_0[1, :l_other])
  else
    dyn.free_surp_boy[τ] =
      freesurp(dyn,
               proj.val[τ-1, :invest],
               proj.val[τ-1, :tpg] + proj.val[τ-1, :l_other])
  end
end

## cashflow projection  -----------------------------------------
"""
`val0!(cap_mkt::CapMkt, invs::InvPort, liabs::LiabIns,
  l_other::LiabOther, proj::Projection)`

Valuation at time `t_0`

**Changed:** `proj::Projection`
"""
function val0!(cap_mkt::CapMkt,
               invs::InvPort,
               liabs::LiabIns,
               l_other::LiabOther,
               proj::Projection)
  proj.val_0[1, :invest] = invs.mv_0
  for mp in liabs.mps
    if 0 <= mp.dur
      proj.val_0[1, :tpg] += tpg(0, cap_mkt.rfr.x, mp)
    end
  end
  proj.val_0[1, :l_other] = pv(0, cap_mkt, l_other)
  proj.val_0[1, :surplus] =
    proj.val_0[1, :invest] -
    proj.val_0[1, :tpg] -
    proj.val_0[1, :l_other]
end

"""
`projectboy!(τ, proj::Projection, liabs::LiabIns)`

Project one year, update values at the beginning of the year `τ`

**Changed:** `proj::Projection`, `liabs::LiabIns  (liabs.mp)`
"""
function projectboy!(τ, proj::Projection, liabs::LiabIns)
  proj.cf[τ, :prem] = 0.0
  proj.cf[τ, :λ_boy] = 0.0
  for mp in liabs.mps
    if τ <= mp.dur
      mp.lx_boy[τ] = (τ == 1 ? 1 : mp.lx_boy_next)
      proj.cf[τ, :prem] += mp.lx_boy[τ] * mp.β[τ, :prem]
      proj.cf[τ, :λ_boy] -=
        mp.lx_boy[τ] * mp.λ[τ, :boy] *
        mp.λ[τ, :cum_infl] / (1 + mp.λ[τ, :infl])
    end
  end
end

"""
`projecteoy!(τ, cap_mkt::CapMkt, invs::InvPort, liabs::LiabIns,
  dyn::Dynamic, proj::Projection)`

Project one year, update values at the end of the year `τ`

**Changed:** `proj::Projection`, `liabs::LiabIns  (liabs.mp)`
"""
function projecteoy!(τ,
                     cap_mkt::CapMkt,
                     invs::InvPort,
                     liabs::LiabIns,
                     dyn::Dynamic,
                     proj::Projection)
  tpg_price_positive = 0.0
  for mp in liabs.mps
    if τ <= mp.dur
      tpg_price_positive +=
        mp.lx_boy[τ] *
        max(0, (τ == 1 ? mp.tpg_price_0 : mp.tpg_price[τ-1]))
    end
  end
  for mp in liabs.mps
    if τ <= mp.dur
      prob = deepcopy(mp.prob)
      prob[:,:sx] *=
        δsx(τ, cap_mkt, invs, mp, dyn)
      prob[:,:px] = 1 .- prob[:,:qx] - prob[:,:sx]
      mp.lx_boy_next = mp.lx_boy[τ] * prob[τ, :px]
      for wx in [:qx, :sx, :px]
        proj.cf[τ, wx] -=
          mp.lx_boy[τ] * prob[τ, wx] * mp.β[τ, wx]
      end
      proj.cf[τ, :λ_eoy] -=
        mp.lx_boy[τ] * mp.λ[τ, :eoy] * mp.λ[τ, :cum_infl]
      proj.val[τ, :tpg] +=
        mp.lx_boy[τ] * prob[τ, :px] *
        tpg(τ, cap_mkt.rfr.x, prob, mp.β, mp.λ)
    end
  end
  proj.cf[τ, :Δtpg] =
    -(proj.val[τ, :tpg] -
        (τ == 1 ? proj.val_0[1, :tpg] : proj.val[τ - 1, :tpg]))

end

"""
`project!(τ, cap_mkt::CapMkt, invs::InvPort, dyn::Dynamic,
  proj::Projection)`

Project one year, investment results from the year `τ`

**Changed:** `proj::Projection`, `invs::InvPort`
"""
function project!(τ,
                  cap_mkt::CapMkt,
                  invs::InvPort,
                  dyn::Dynamic,
                  proj::Projection)
  mv_boy =
    (τ == 1 ? proj.val_0[1, :invest] : proj.val[τ - 1, :invest])
  mv_boy += proj.cf[τ, :prem] + proj.cf[τ, :λ_boy]
  alloc!(τ, cap_mkt, invs)
  project!(τ, mv_boy, invs)
  proj.cf[τ, :invest] = invs.mv[τ] - mv_boy
end

"""
`bonus!(τ, invs::InvPort, liabs::LiabIns, dyn::Dynamic,
  proj, surp_pre_profit_tax_bonus)`

Bonus at the end of year `τ`

**Changed:** `proj::Projection`
"""
function bonus!(τ,
                invs::InvPort,
                liabs::LiabIns,
                dyn::Dynamic,
                proj,
                surp_pre_profit_tax_bonus)
  for mp in liabs.mps
    if τ <= mp.dur
      proj.cf[τ, :bonus] -=
        min(surp_pre_profit_tax_bonus,
            mp.lx_boy[τ] *
              bonusrate(τ, invs.yield[τ], mp, dyn) *
              max(0, (τ == 1 ?
                        mp.tpg_price_0 :
                        mp.tpg_price[τ-1])))
    end
  end
end

"""
`investpredivid(τ, invs::InvPort, proj::Projection)`

Market value of assets before payment of dividends
"""
function investpredivid(τ, invs::InvPort, proj::Projection)
  invs.mv_boy[τ] +
    sum(convert(Array,
                proj.cf[τ,
                        [:invest, :qx, :sx, :px, :λ_eoy, :bonus,
                         :l_other, :tax, :gc]]))
end

"""
`project!(τ, cap_mkt::CapMkt, invs::InvPort, liabs::LiabIns,
  liab_other::LiabOther, dyn::Dynamic, proj::Projection)`

Project one year

**Changed:** `proj::Projection`, `invs::InvPort`, `dyn::Dynamic`
"""
function project!(τ,
                  cap_mkt::CapMkt,
                  invs::InvPort,
                  liabs::LiabIns,
                  liab_other::LiabOther,
                  dyn::Dynamic,
                  proj::Projection)
  projectboy!(τ, proj, liabs)
  proj.cf[τ, :new_debt] = -getloan(τ, liab_other)
  project!(τ, cap_mkt, invs, dyn, proj)
  update!(τ, proj, dyn)
  proj.cf[τ, :λ_eoy] = -invs.cost[τ]
  proj.cf[τ, :l_other] = -paycoupon(τ, liab_other)
  proj.cf[τ, :l_other] -= payprincipal(τ, liab_other)
  proj.val[τ, :l_other] = pv(τ, cap_mkt, liab_other)
  projecteoy!(τ, cap_mkt, invs, liabs, dyn, proj)

  proj.cf[τ,:tax] = 0.0
  proj.cf[τ,:bonus] = 0.0
  surp_pre_profit_tax_bonus =
    max(0,
        investpredivid(τ, invs, proj) -
          proj.val[τ, :tpg] -
          proj.val[τ, :l_other]  )
  bonus!(τ, invs, liabs, dyn, proj, surp_pre_profit_tax_bonus)
  proj.cf[τ, :profit] =
    sum(convert(Array, proj.cf[τ, [:prem, :invest,
                                   :qx, :sx, :px, :λ_boy, :λ_eoy,
                                   :Δtpg, :bonus, :l_other]]))
  tax = proj.tax_rate * proj.cf[τ, :profit] ## could be negative
  tax_credit_pre =
    (τ == 1 ? proj.tax_credit_0 : proj.tax_credit[τ - 1])
  proj.cf[τ, :tax] = -max(0, tax - tax_credit_pre)
  proj.tax_credit[τ] =  ## no new tax credit generated
  tax_credit_pre - tax - proj.cf[τ, :tax]

  proj.val[τ, :invest] = investpredivid(τ, invs, proj)
  proj.cf[τ, :divid] =
    -freesurp(dyn,
              proj.val[τ, :invest],
              proj.val[τ, :tpg] + proj.val[τ, :l_other])
  proj.val[τ, :invest] +=  proj.cf[τ, :divid]

  proj.val[τ, :surplus] =
    proj.val[τ, :invest] -
    proj.val[τ, :tpg] -
    proj.val[τ, :l_other]
end

"""
`pvprev(rfr, cf, pv) = (cf + pv) /  (1 + rfr)`

Recursive step in generic present value calculation
"""
pvprev(rfr, cf, pv) = (cf + pv) /  (1 + rfr)

"""
`pvvec(rfr::Vector{Float64}, cf)``

Generic present value calculation, where `rfr` denotes the risk
free rate for each year and `cf`. The length of the vector `rfr`
may not be smaller than the length of the cashflow.
Each payment occurs at the end of the corresponding year.

Output: A vector of the same length as `cf` contaning the present
value at the end of each year. The last component is zero.
"""
function pvvec(rfr::Vector{Float64}, cf)
  T = length(cf)
  val = zeros(Float64, T)
  for t in reverse(collect(1:(T-1))) # [T-1:-1:1] #r
    val[t] = pvprev(rfr[t + 1], val[t + 1], cf[t + 1])
  end
  return val
end

"""
`valbonus!(rfr::Vector{Float64}, proj::Projection) `

Provisions for future bonus payments (at each time `τ`).
Needs to be called after the projection is completed

**Changed:** proj::Projection
"""
function valbonus!(rfr::Vector{Float64},
                   proj::Projection)
  proj.val[:bonus] = pvvec(rfr,  -proj.cf[:bonus])
  proj.val_0[:bonus] =
    pvprev(rfr[1], -proj.cf[1, :bonus], proj.val[1, :bonus])
end

"""
`valcostprov!(rfr::Vector{Float64}, invs::InvPort,
  proj::Projection)`

Provisions for future costs (at each time `τ`). Absolute costs
(including absolute investment costs from *all* investments) and
relative investment costs for provisions are considered.
It is assumed that provisions are backed by cash investments
Needs to be called after the projection is completed

**Changed:** proj::Projection
"""
function valcostprov!(rfr::Vector{Float64},
                      invs::InvPort,
                      proj::Projection)
  cash_cost = deepcopy(invs.igs[:IGCash].cost)
  proj.cf[1, :cost_prov] =
    proj.fixed_cost_gc[1] +
    sum(convert(Array,
                proj.val_0[1, [:tpg, :bonus, :l_other]])) *
    cash_cost.cum_infl_rel[1] * cash_cost.rel[1]
  for t = 2:proj.dur
    proj.cf[t, :cost_prov] =
      proj.fixed_cost_gc[t] +
      sum(convert(Array,
                  proj.val[t - 1, [:tpg, :bonus, :l_other]])) *
      cash_cost.cum_infl_rel[t] * cash_cost.rel[t]
  end
  proj.val[:cost_prov] =
    pvvec(rfr - cash_cost.rel,  proj.cf[:cost_prov])
  proj.val_0[:cost_prov] = pvprev(rfr[1] - cash_cost.rel[1],
                                  proj.cf[1, :cost_prov],
                                  proj.val[1, :cost_prov])
end
