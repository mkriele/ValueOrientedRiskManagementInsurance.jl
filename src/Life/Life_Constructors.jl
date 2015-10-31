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

function ModelPoint(n, t_0, t_start,
                    prob_be, sx_be_fac, λ_be,
                    cost_infl,
                    hypo_bonus_rate, product, ins_sum,
                    pension_contract)
  ## Time model for model point: See documentation of type
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
                    ones(Float64, dur),
                    pension_contract)
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
                          df_port[d, :ins_sum],
                          df_port[d, :pension_contract]))
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
  valcostprov!(cap_mkt.rfr.x, invs, proj)
  return proj
end

