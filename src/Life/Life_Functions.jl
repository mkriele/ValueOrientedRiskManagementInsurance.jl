export eqshock
export goingconcern!

## constructors =================================================

## capital market -----------------------------------------------
CapMkt(proc_stock, proc_rfr) =
  CapMkt(deepcopy(proc_stock), deepcopy(proc_rfr))
## assets -------------------------------------------------------
function IGCost(df_cost)
  IGCost(df_cost[:, :rel],
         df_cost[:, :abs],
         df_cost[:, :infl_rel],
         df_cost[:, :infl_abs],
         cumprod(1 .+ df_cost[:, :infl_rel]),
         cumprod(1 .+ df_cost[:, :infl_abs]),
         zeros(Float64, nrow(df_cost))
         )
end

function IGStock(cap_mkt::CapMkt, mv_0, alloc, cost)
  investments = Array(InvestStock, 0)
  for i = 1:length(alloc.name)
    push!(investments,
          InvestStock(alloc.name[i],
                      cap_mkt.stock,
                      mv_0 * alloc.total[1] * alloc.all[1,i],
                      zeros(Float64, size(alloc.all, 1))))
  end
  IGStock(investments, mv_0 * alloc.total[1],
          zeros(Float64, size(alloc.all, 1)), alloc,
          deepcopy(cost))
end

function IGCash(cap_mkt::CapMkt, mv_0, alloc, cost)
  investments = Array(InvestCash, 0)
  for i = 1:length(alloc.name)
    push!(investments,
          InvestCash(alloc.name[i],
                     cap_mkt.rfr,
                     mv_0 * alloc.total[1] * alloc.all[1,i],
                     zeros(Float64, size(alloc.all, 1)),
                     alloc.lgd[i],
                     alloc.cqs[i]))
  end
  IGCash(investments, mv_0 * alloc.total[1],
         zeros(Float64, size(alloc.all, 1)), alloc,
         deepcopy(cost))
end

function InvPort(t_0,
                 dur,
                 cap_mkt::CapMkt,
                 mv_0,
                 allocs::Dict{Symbol, Alloc},
                 costs::Dict{Symbol, DataFrame}
                 )
  igs = Dict{Symbol, InvestGroup}()
  for ig_symb in collect(keys(allocs))
    ## ig_symb are the symbols corresponding to the
    ## types of investment groups: :IGCash, IGStock
    merge!(igs, [ig_symb => eval(ig_symb)(cap_mkt,
                                          mv_0,
                                          allocs[ig_symb],
                                          IGCost(costs[ig_symb])
                                          )])
  end
  return(InvPort(t_0,                  ## start of projection
                 mv_0,                 ## init. mv pre prem..
                 zeros(Float64, dur),  ## mv bop post prem...
                 zeros(Float64, dur),  ## mv eop
                 zeros(Float64, dur),  ## average yield
                 zeros(Float64, dur),  ## investment costs
                 igs))                 ## investment groups
end

## insurance liabilities ----------------------------------------
function Product(rfr_price, prob_price, β_in, λ_price)
  dur = nrow(β_in)
  prob = deepcopy(prob_price)
  prob[:px] = 1 .- prob[:qx] - prob[:sx]
  λ_price[:cum_infl] = cumprod(1 .+ λ_price[:infl])
  return Product(dur, rfr_price, prob, β_in, λ_price,
                 premium(1, rfr_price, prob, β_in, λ_price))
end

## Time model for model point:
##
## project time τ:            0                   dur
## real time t:      t_start  t_0
## product time s:   0        s_0                 product.dur
##                   |--------|-------------------|--------------
##                            \-------------------/
##                                     dur
##                   \----------------------------/
##                           product.dur
function ModelPoint(n, t_0, t_start,
                    prob_be, sx_be_fac, λ_be,
                    cost_infl,
                    hypo_bonus_rate, product, ins_sum)
  s_0 = t_0 - t_start
  dur = product.dur - s_0
  s_future = (s_0 + 1):product.dur
  prob = deepcopy(prob_be)[s_future, :]
  prob[:sx] *= sx_be_fac
  prob[:px] = 1 .- prob[:qx] - prob[:sx]
  lx_boy = zeros(Float64, dur)
  β = DataFrame()
  for name in names(product.β)
    β[name] = n * ins_sum * product.β[s_future, name]
  end
  β[:prem] .*= product.prem_norm
  β[:sx] .*= product.prem_norm

  λ = deepcopy(λ_be)[s_future, :]
  λ[:boy] .*= n * ins_sum
  λ[:eoy] .*= n * ins_sum
  ## be cost inflation input relates to t_0 not s_0:
  λ[:infl] = deepcopy(cost_infl)
  λ[:cum_infl] = cumprod(1 .+ λ[:infl])
  λ_price = deepcopy(product.λ)[s_future, :]
  λ_price[:boy] .*= n * ins_sum
  λ_price[:eoy] .*= n * ins_sum
  λ_price[:cum_infl] = cumprod(1 .+ product.λ[:infl])[s_future]
  rfr_price_0 = product.rfr[s_0 == 0 ? 1 : s_0]
  rfr_price = product.rfr[s_future]
  tpg_price_0 = tpg(0,
                    rfr_price,
                    product.prob[s_future, :],
                    β,
                    λ_price)
  tpg_price = zeros(Float64, dur)
  for τ = 1:dur
    tpg_price[τ] = tpg(τ,
                       rfr_price,
                       product.prob[s_future, :],
                       β,
                       λ_price)
  end
  return ModelPoint(n, t_start, dur, prob, lx_boy, 0.0,
                    β, λ,  hypo_bonus_rate,
                    rfr_price_0, rfr_price,
                    tpg_price_0, tpg_price,
                    ones(Float64, dur))
end

function LiabIns(t_0, prob_be, λ_be, cost_infl, product, df_port)
  n = nrow(df_port)
  mps = Array(ModelPoint, 0)
  dur = 0
  for d = 1:n
    push!(mps, ModelPoint(df_port[d, :n],
                          t_0,
                          df_port[d, :t_start],
                          prob_be,
                          df_port[d, :sx_be_fac],
                          λ_be,
                          cost_infl[d],
                          df_port[d, :bonus_rate_hypo],
                          product,
                          df_port[d, :ins_sum]))
    dur = max(dur, mps[d].dur)
  end
  gc = zeros(Float64, dur)
  for mp in mps
    mp.gc = zeros(Float64, dur)
    mp.gc[1:mp.dur] +=
      vcat(1, cumprod(mp.prob[1:(mp.dur-1), :px]))
    gc +=  mp.n * mp.gc
  end
  gc ./= gc[1]
  Δgc = diff(vcat(gc, 0))
  return LiabIns(n, t_0, dur, mps, gc, Δgc)
end

## other liabilities --------------------------------------------
function Debt(t_0, t_debt_0, t_debt_mat,
              name::Symbol, nominal, coupon)
  name = name
  τ_debt_0 = t_debt_0 - t_0
  dur_debt = t_debt_mat - t_debt_0 + 1
  τ_mat = dur_debt + τ_debt_0 - 1
  Debt(name, t_debt_0, t_debt_mat, τ_debt_0, τ_mat, nominal, coupon)
end

function Debt(t_0, df_debt::DataFrame)
  name = df_debt[1, :name]
  t_init = df_debt[1, :t_init]
  t_mat = df_debt[1, :t_mat]
  τ_init = t_init - t_0
  dur = t_mat - t_init + 1
  τ_mat = dur + τ_init - 1
  Debt(name, t_init, t_mat, τ_init, τ_mat,
       df_debt[1, :nominal], df_debt[1, :coupon])
end

function LiabOther(t_0, df_debts::DataFrame)
  subord = Array(Debt, nrow(df_debts))
  for d = 1:nrow(df_debts)
    subord[d] = Debt(t_0, df_debts[d,:])
  end
  return LiabOther(subord)
end

## dynamics -----------------------------------------------------
function Dynamic(dur, bonus_factor, quota_surp)
  Dynamic(bonus_factor,
          quota_surp,
          zeros(Float64, dur))
end

## cashflow projection ------------------------------------------
function Projection(liabs, tax_rate, tax_credit_0)
  t_0 = liabs.t_0
  dur = liabs.dur
  cf = DataFrame(
    qx = zeros(Float64, dur),
    sx = zeros(Float64, dur),
    px = zeros(Float64, dur),
    prem = zeros(Float64, dur),
    λ_boy = zeros(Float64, dur),
    λ_eoy = zeros(Float64, dur),
    bonus = zeros(Float64, dur),
    invest = zeros(Float64, dur),
    new_debt = zeros(Float64, dur),
    l_other = zeros(Float64, dur),
    profit = zeros(Float64, dur),
    tax = zeros(Float64, dur),
    profit = zeros(Float64, dur),
    divid = zeros(Float64, dur),
    gc = zeros(Float64, dur),        ## not affecting profit/loss
    Δtpg = zeros(Float64, dur),      ## not real cf, affects p/l
    cost_prov = zeros(Float64, dur)  ## for the cost provisions
    )
  val = DataFrame(
    invest = zeros(Float64, dur),
    tpg = zeros(Float64, dur),
    l_other = zeros(Float64, dur),
    surplus = zeros(Float64, dur),
    bonus = zeros(Float64, dur),
    cost_prov = zeros(Float64, dur)      ## cost provisions
    )
  val_0 = deepcopy(val[1, :])
  return Projection(t_0, dur, cf, val_0, val,
                    tax_rate,
                    tax_credit_0,
                    zeros(Float64, dur),
                    zeros(Float64, dur))
end

function Projection(tax_rate,
                    tax_credit_0,
                    cap_mkt::CapMkt,
                    invs::InvPort,
                    liabs::LiabIns,
                    liabs_other::LiabOther,
                    dyn::Dynamic)
  proj = Projection(liabs, tax_rate, tax_credit_0)
  for τ = 1:liabs.dur
    for ig in [:IGCash, :IGStock]
      proj.fixed_cost_gc[τ] +=
        invs.igs[ig].cost.abs[τ] *
        invs.igs[ig].cost.cum_infl_abs[τ] *
        liabs.gc[τ]
    end
  end
  l_other = deepcopy(liabs_other)
  goingconcern!(l_other, liabs.Δgc)
  val0!(cap_mkt, invs, liabs, l_other, proj)
  proj.cf[:,:gc] = liabs.Δgc * (proj.val_0[1,:invest] -
                                  proj.val_0[1,:tpg] -
                                  proj.val_0[1,:l_other])
  for τ = 1:liabs.dur
    project!(τ, cap_mkt, invs, liabs, l_other, dyn, proj)
  end
  valbonus!(cap_mkt.rfr.x, proj)
  valcostprov!(cap_mkt.rfr.x, invs.igs[:IGCash].cost, proj)
  return proj
end

## other functions ##############################################
## capital market -----------------------------------------------
function yield(τ, stock::Stock)
  if τ > 1
    return stock.x[τ] / stock.x[τ - 1] - 1
  elseif τ == 1
    return stock.x[τ] / stock.x_0 - 1
  else
    return stock.yield_0
  end
end

function yield(τ, rfr::RiskFreeRate)
  if τ >= 1
    return rfr.x[τ]
  else
    return rfr.yield_0
  end
end

yieldmkt(τ, cap_mkt::CapMkt) = yield(τ, cap_mkt.stock)

yieldrfr(τ, cap_mkt::CapMkt) = yield(τ, cap_mkt.rfr)

## (1+f[1])(1+f[2])...(1+f[n]) = (1+s[n])^n
forw2spot(f::Vector{Float64}) =
  cumprod(1 .+ f) .^ (1 ./ [1:length(f)]) -1

## (1 + s[n-1])^(n-1) * (1+f[n]) = (1+s[n])^n
function spot2forw(s::Vector{Float64})
  f = zeros(Float64, length(s))
  for n in length(s):-1:2
    f[n] = (1 + s[n])^n / (1 + s[n-1])^(n - 1) -1
  end
  f[1] = s[1]
  return f
end

## investments --------------------------------------------------
function project!(τ::Int,
                  mv_bop::Float64,
                  invest::Invest)            ## changed
  invest.mv[τ] = (1 + yield(τ, invest.proc)) * mv_bop
end

function project!(τ::Int,
                  mv_bop_total::Float64,  ## total over all igs
                  ig::InvestGroup)        ## changed
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

function project!(τ::Int,
                  cap_mkt::CapMkt,
                  mv_boy::Float64,
                  invs::InvPort)          ## changed
  alloc!(τ, cap_mkt, invs)
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
function premium(ins_sum, rfr, prob, β, λ)
  lx_boy = [1,cumprod(prob[:px])[1:end-1]]
  v_eoy = 1 ./ cumprod(1 .+ rfr)
  v_boy = [1, v_eoy[1:end-1]]
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

## technical provisions (recursion)
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

tpg(τ, rfr, mp) =
  tpg(τ, rfr, mp.prob, mp.β, mp.λ)


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

function pv(τ::Int, cap_mkt::CapMkt, l_other::LiabOther)
  p_v = 0.0
  for debt in l_other.subord
    p_v += pv(τ, cap_mkt, debt)
  end
  return p_v
end

paycoupon(τ::Int, debt::Debt) =
  (debt.τ_init <= τ <= debt.τ_mat ? debt.coupon : 0.0)

function paycoupon(τ::Int, l_other::LiabOther)
  pay = 0.0
  for debt in l_other.subord
    pay += paycoupon(τ, debt)
  end
  return pay
end

payprincipal(τ::Int, debt::Debt) =
  (τ == debt.τ_mat ? debt.nominal : 0.0)

function payprincipal(τ::Int, l_other::LiabOther)
  pay = 0.0
  for debt in l_other.subord
    pay += payprincipal(τ, debt)
  end
  return pay
end

getloan(τ::Int, debt::Debt) =
  (τ == debt.τ_init ? debt.nominal : 0.0)

function getloan(τ::Int, l_other::LiabOther)
  nominal = 0.0
  for debt in l_other.subord
    nominal += getloan(τ, debt)
  end
  return nominal
end

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

function goingconcern!(l_other::LiabOther, Δgc::Vector{Float64})
  l_other.subord = goingconcern(l_other.subord, Δgc)
end


## dynamics -----------------------------------------------------
## state of the economy
dynstate(τ, cap_mkt::CapMkt) =
  yieldmkt(τ, cap_mkt) / max(yieldrfr(τ, cap_mkt), eps()) -1

## averaged state of the economy (two years, only for stocks)
dynstateavg(τ, cap_mkt::CapMkt) =
  0.5 * (yieldmkt(τ - 1, cap_mkt) /
           max(yieldrfr(τ-1, cap_mkt), eps())
         + yieldmkt(τ, cap_mkt) /
           max(yieldrfr(τ, cap_mkt), eps())) - 1

## dynamic update of asset allocations
function alloc!(τ,
                cap_mkt::CapMkt,
                invs::InvPort)             ## changed
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

## dynamic bonus rate declaration
bonusrate(τ, yield_eoy, rfr_price, bonus_factor) =
  max(bonus_factor * (yield_eoy - rfr_price), 0.0)

bonusrate(τ, yield_eoy, mp::ModelPoint, dyn) =
  bonusrate(τ,
            yield_eoy,
            ( τ == 0 ? mp.rfr_price_0 : mp.rfr_price[τ]),
            dyn.bonus_factor)

function biquotient(τ, yield_eoy, cap_mkt,
                    invs::InvPort, mp, dyn)
  ## Compare currend yield with bonus expectation
  ## based on experience from previous year
  ind_bonus =
    yieldmkt(τ, cap_mkt) /
    max(eps(),
        bonusrate(τ - 1, yield_eoy, mp, dyn) +
          yieldrfr(τ, cap_mkt))
  ind_bonus_hypo =
    yieldmkt(0, cap_mkt) /
    max(eps(), mp.bonus_rate_hypo + yieldrfr(0, cap_mkt))
  return ind_bonus / ind_bonus_hypo
end

## dynamic lapse probabilities (factor to adjust init. estimate)
function δsx(τ, cap_mkt::CapMkt, invs::InvPort, mp::ModelPoint,
             dyn)
  yield_eoy =
    invs.igs[:IGStock].alloc.total[τ] * yieldmkt(τ, cap_mkt) +
    invs.igs[:IGCash].alloc.total[τ] * yieldrfr(τ, cap_mkt)
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

## dynamic dividend declaration
freesurp(dyn, invest_pre, liab) =
 max(0, invest_pre - (1 + dyn.quota_surp) * liab)

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


## update dynamic parameters
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
function val0!(cap_mkt::CapMkt,
               invs::InvPort,
               liabs::LiabIns,
               l_other::LiabOther,
               proj::Projection)       ## changed
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

function project!(τ,
                  cap_mkt::CapMkt,
                  invs::InvPort,     ## changed
                  dyn::Dynamic,  ## changed
                  proj::Projection)      ## changed
  mv_boy =
    (τ == 1 ? proj.val_0[1, :invest] : proj.val[τ - 1, :invest])
  mv_boy += proj.cf[τ, :prem] + proj.cf[τ, :λ_boy]
  project!(τ, cap_mkt, mv_boy, invs)
  proj.cf[τ, :invest] = invs.mv[τ] - mv_boy
end

function bonus!(τ,
                invs::InvPort,
                liabs::LiabIns,
                dyn::Dynamic,
                proj,              ## changed
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

function investpredivid(τ,
                        invs::InvPort,
                        proj::Projection)
  invs.mv_boy[τ] +
    sum(array(proj.cf[τ,
                      [:invest, :qx, :sx, :px, :λ_eoy, :bonus,
                       :l_other, :tax, :gc]]))
end

function project!(τ,
                  cap_mkt::CapMkt,
                  invs::InvPort,     ## changed
                  liabs::LiabIns,
                  liab_other::LiabOther,
                  dyn::Dynamic,      ## changed
                  proj::Projection  ## changed
                  )
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
    sum(array(proj.cf[τ, [:prem, :invest,
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

pvprev(rfr, cf, pv) = (cf + pv) /  (1 + rfr)

function pvvec(rfr::Vector{Float64}, cf)
  T = length(cf)
  val = zeros(Float64, T)
  for t in [T-1:-1:1]
    val[t] = pvprev(rfr[t + 1], val[t + 1], cf[t + 1])
  end
  return val
end

function valbonus!(rfr::Vector{Float64},
                   proj::Projection)        ## changed
  proj.val[:bonus] = pvvec(rfr,  -proj.cf[:bonus])
  proj.val_0[:bonus] =
    pvprev(rfr[1], -proj.cf[1, :bonus], proj.val[1, :bonus])
end

function valcostprov!(rfr::Vector{Float64},
                      cost::IGCost,
                      proj::Projection)       ## changed
  proj.cf[1, :cost_prov] =
    proj.fixed_cost_gc[1] +
    sum(array(proj.val_0[1, [:tpg, :bonus, :l_other]])) *
    cost.cum_infl_rel[1] * cost.rel[1]
  for t = 2:proj.dur
    proj.cf[t, :cost_prov] =
      proj.fixed_cost_gc[t] +
      sum(array(proj.val[t - 1, [:tpg, :bonus, :l_other]])) *
      cost.cum_infl_rel[t] * cost.rel[t]
  end
  proj.val[:cost_prov] =
    pvvec(rfr - cost.rel,  proj.cf[:cost_prov])
  proj.val_0[:cost_prov] = pvprev(rfr[1] - cost.rel[1],
                                  proj.cf[1, :cost_prov],
                                  proj.val[1, :cost_prov])
end
