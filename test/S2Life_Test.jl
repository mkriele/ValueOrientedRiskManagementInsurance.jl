using ValueOrientedRiskManagementInsurance
using DataFrames
using Base.Test

include("S2Life.jl")

VORMI = ValueOrientedRiskManagementInsurance

println("Start S2Life test")

ins_sum = df_portfolio[1, :ins_sum]

for t = 1: T
  @test prob_price[1, :sx] == prob_price[t, :sx]
  @test rfr_price[1] == rfr_price[t]
  @test λ_price[1, :infl] == λ_price[t, :infl]
  @test λ_price[1, :eoy] == λ_price[t, :eoy]
  if t > 1
    @test λ_price[t, :boy] == 0
  end
  @test  convert(Array, β[:sx]) == cumsum(fill(β[1, :sx], 5))
  ## following is used in presentation of pricing calculation:
  @test prob_price[t, :qx] ≈ (10 + t -1)/10000
end

for i = 1:nrow(df_portfolio)
  @test df_portfolio[i, :ins_sum] == ins_sum
end


## Premium ------------------------------------------------------
prob_price[:px] = 1 .- prob_price[:, :qx] - prob_price[1, :sx]
lx_price_boy =
  convert(Array, cumprod(prob_price[:px]) ./ prob_price[:px])

v_price_eoy = cumprod(1 ./ (1 .+ rfr_price))
v_price_boy = v_price_eoy .* (1 .+ rfr_price)

infl_price_eoy = convert(Array, cumprod(1+λ_price[:infl]))
infl_price_boy =
  infl_price_eoy ./ convert(Array, 1+λ_price[:infl])

prem_price_ratio =
  sum(lx_price_boy .* v_price_boy .*
      λ_price[:, :boy] .* infl_price_boy +
        lx_price_boy .* v_price_eoy .*
      (λ_price[:, :eoy] .* infl_price_eoy +
         prob_price[:, :qx] .* β[:, :qx] +
         prob_price[:, :px] .* β[:, :px])) /
  sum(lx_price_boy .* v_price_boy -
        prob_price[:, :sx] .*
      β[:, :sx] .* lx_price_boy .* v_price_eoy)

for i = 1:nrow(df_portfolio)
  for t = 1:i
    @test prem_price_ratio ≈
                    liab_ins.mps[i].β[t, :prem] ./liab_ins.mps[i].β[t, :qx]
  end
end

prem_price = prem_price_ratio * ins_sum

## techn. prov. (pricing) calc ----------------------------------
tp_price = zeros(Float64, T)
for t = (T-1):-1:1
  tp_price[t] =
    - prem_price +
    1 / (1 + rfr_price[t + 1]) *
    (λ_price[t + 1, :eoy] * infl_price_eoy[t + 1] *ins_sum +
       prob_price[t+1, :qx] * ins_sum +
       prob_price[t + 1, :sx] * β[t + 1, :sx] * prem_price +
       prob_price[t + 1, :px] * (β[t + 1, :px] * ins_sum +
                                   tp_price[t + 1]))
end

for i = 1:nrow(df_portfolio)
  t_contract = T + df_portfolio[i, :t_start]
  for t = 1:t_contract
    @test tp_price[T - t_contract + 1 : T] ≈
          liab_ins.mps[i].tpg_price / df_portfolio[i,:n]
  end
end

## best estimate assumptions ------------------------------------

y_stock =
  vcat(proc_stock.x[1]/proc_stock.x_0 - 1,
       Float64[proc_stock.x[t] / proc_stock.x[t-1] - 1 for t in 2:T])
delta_qx = prob_price[1, :qx] - prob_be[1, :qx]
for t = 1:T
  @test y_stock[t] ≈ proc_stock.x[1]-1
end
for t = 1:T
  @test prob_be[t, :qx] + delta_qx == prob_price[t, :qx]
end

## Costs ========================================================

## In our example the cost inflation is constant in time
for d = 1:nrow(df_portfolio)
  for x in cost_infl_be[d]
    @test x == cost_infl_be[1][1]
  end
end
## In our example all costs are constant in time
for t = 1:T
  if t > 1
    @test  λ_be[t, :boy] == 0
  end
  @test λ_be[1, :eoy] == λ_be[t, :eoy]
  @test λ_invest[:IGCash][1, :rel] == λ_invest[:IGCash][t, :rel]
  @test λ_invest[:IGCash][1, :abs] == λ_invest[:IGCash][t, :abs]
  @test λ_invest[:IGStock][1, :rel] == λ_invest[:IGStock][t, :rel]
  @test λ_invest[:IGStock][1, :abs] == λ_invest[:IGStock][t, :abs]
end

## State of the economy =========================================

## notice that indices of tmp_stock, tmp_state are shifted by one
tmp_stock = [cap_mkt.stock.x_0; cap_mkt.stock.x]
state = Float64[(tmp_stock[t+1]/tmp_stock[t]-1-rfr[t])/rfr[t]
                for t = 1:T]
tmp_state =
  [cap_mkt.stock.yield_0 / cap_mkt.rfr.yield_0 - 1; state]
state_avg =
  Float64[(tmp_state[t + 1] + tmp_state[t]) / 2 for t = 1:T]


allocation = 0.5 * (1 .- exp.(-max.(0, state_avg)))
## restore initial allocation
allocation[1] = invs.igs[:IGStock].alloc.total[1]

state_orig =
  Float64[VORMI.dynstate(t,cap_mkt) for t in 1:T]
state_avg_orig =
  Float64[VORMI.dynstateavg(t,cap_mkt) for t in 1:T]
@test state ≈ state_orig
@test state_avg ≈ state_avg_orig
@test allocation ≈ invs.igs[:IGStock].alloc.total

y_invest =
  Float64[allocation[t] * y_stock[t] + (1-allocation[t]) *rfr[t]
          for t in 1:T]

t_bonus_quota = dyn.bonus_factor * (y_invest - rfr_price)

sx_basis = Array{Vector{Float64}}(nrow(df_portfolio))
for i = 1:nrow(df_portfolio)
  sx_basis[i] = convert(Array, liab_ins.mps[i].prob[:sx])
end

v_eoy = cumprod(1 ./ (1 .+ rfr))
v_boy = v_eoy .* (1 .+ rfr)

infl_eoy = cumprod(1 + cost_infl_be[end])
infl_boy = infl_eoy ./ (1 .+ cost_infl_be[end])

prob = deepcopy(prob_be)
prob[:px] = 1 .- prob[:qx] - prob[:sx]

rfr_cost = rfr - λ_invest[:IGCash][:, :rel]

function tpberec(tp_next, t, τ, prob_sx)
  prob_px = 1.- prob[:qx] - prob_sx
  - prem_price +
    λ_be[t + τ - t_0 + 1, :boy] * infl_boy[t + 1] * ins_sum+
    1 / (1 + rfr[t + 1]) *
    (λ_be[t + τ - t_0 + 1, :eoy] * infl_eoy[t + 1] * ins_sum +
       prob[t + τ - t_0 + 1, :qx] * ins_sum +
       prob_sx[t + τ - t_0 + 1] * β[t + τ - t_0 + 1, :sx] * prem_price +
       prob_px[t + τ - t_0 + 1] * (β[t + τ - t_0 + 1, :px] * ins_sum + tp_next))
end

τ=1
d=4
t=3
prem_price
mp = deepcopy(liab_ins.mps[d])
fn = df_portfolio[d, :n]
prob_sx = convert(Array, prob[:sx]) * df_portfolio[d, :sx_be_fac]
prob_px = 1.- prob[:qx] - prob_sx
@test fn * prem_price ≈ mp.β[t+1, :prem]
@test mp.λ[t + 1, :boy] *  mp.λ[t + 1, :cum_infl] /
      (1 + mp.λ[t + 1, :infl]) ≈
      λ_be[t + τ - t_0 + 1, :boy] * infl_eoy[t + 1] * ins_sum
disc = 1/(1 + rfr[t + 1])
@test 1 / (1 + rfr_cost[t + 1]) ≈
      disc/(1 - disc * invs.igs[:IGCash].cost.rel[τ + 1])
@test fn * λ_be[t + τ - t_0 + 1, :eoy] *
      infl_eoy[t + 1] * ins_sum ≈
      mp.λ[t + 1, :eoy] * mp.λ[t + 1, :cum_infl]
@test fn * λ_be[t + τ - t_0 + 1, :eoy] * ins_sum ≈
      mp.λ[t + 1, :eoy]
@test fn * prob[t + τ - t_0 + 1, :qx] * ins_sum ≈
      mp.prob[t + 1, :qx] * mp.β[t + 1, :qx]
@test fn * prob_sx[t + τ - t_0 + 1] *
      β[t + τ - t_0 + 1, :sx] * prem_price  ≈
      mp.prob[t + 1, :sx] * mp.β[t + 1, :sx]
@test fn * prob_px[t + τ - t_0 + 1] *
      β[t + τ - t_0 + 1, :px] * ins_sum ≈
      mp.prob[t + 1, :px] * (mp.β[t + 1, :px] + 0)


## Going concern ================================================
tmp_gc = zeros(Float64, T)
for d = 1:nrow(df_portfolio)
  lx_boy = df_portfolio[d, :n]
  tmp_gc[1] += lx_boy
  for t = 2:nrow(liab_ins.mps[d].prob)
    lx_boy *=
      (1 - liab_ins.mps[d].prob[t - 1,:qx] -
         liab_ins.mps[d].prob[t - 1,:sx])
    tmp_gc[t] += lx_boy

  end
end
tmp_gc /= sum(df_portfolio[:n])

@test tmp_gc ≈ liab_ins.gc

tmp_gc_extension = [liab_ins.gc; 0]
tmp_Δgc = zeros(Float64,T)
for t = 1:T
  tmp_Δgc[t] = tmp_gc_extension[t+1] - tmp_gc_extension[t]
end
@test liab_ins.Δgc ≈ tmp_Δgc

## Going concern absolute costs ---------------------------------
cost_abs =
  Float64[liab_ins.gc[t] *
            (λ_invest[:IGCash][t, :abs] *
               prod(1 + (λ_invest[:IGCash][1:t, :infl_abs])) +
               λ_invest[:IGStock][t, :abs] *
               prod(1 + (λ_invest[:IGStock][1:t, :infl_abs])))
          for t = 1:T]

@test proj.fixed_cost_gc ≈ cost_abs

## Subordinated debt --------------------------------------------
# cf_liab_other_unscaled =
#   fill(-df_sub_debt[1, :coupon], df_sub_debt[1, :t_mat] - t_0)
# cf_liab_other_unscaled[df_sub_debt[1, :t_mat] - t_0] -=
#   df_sub_debt[1, :nominal]

l_other = deepcopy(liab_other)
VORMI.goingconcern!(l_other, liab_ins.Δgc)
cf_l_other = Array{Vector{Float64}}(T)
for t = 1:T
  cf_l_other[t] = zeros(Float64, l_other.subord[t].τ_mat)
  fill!(cf_l_other[t], -l_other.subord[t].coupon)
  cf_l_other[t][l_other.subord[t].τ_mat] -=
    l_other.subord[t].nominal
end

cf_l_other_total = zeros(Float64, T)
for t = 1:T
  cf_l_other_total[1:t] += cf_l_other[t]
end

@test sum([cf_l_other[t][end] for t = 1:T]) ≈
      -df_sub_debt[1, :coupon] - df_sub_debt[1, :nominal]
for t = 2:T
  @test cf_l_other[t][1] / cf_l_other[t][end] ≈
        df_sub_debt[1, :coupon] /
        (df_sub_debt[1, :coupon] + df_sub_debt[1, :nominal])
end
@test sum(-[cumprod(1 ./ (1 .+ rfr))[1:t] ⋅ cf_l_other[t]
            for t = 1:T]) ≈
      proj.val_0[1, :l_other]

for d = 1:(T-1)
  @test -cumprod(1 ./ (1 .+ rfr[d+1:T])) ⋅
        cf_l_other_total[d+1:T] ≈
        proj.val[d, :l_other]
end

## gc surplus adjustment ----------------------------------------
@test proj.cf[:gc] ≈
      (proj.val_0[1, :invest] - proj.val_0[1, :tpg] - proj.val_0[1, :l_other]) * liab_ins.Δgc

## Technical provisions for guaranteed benefits, t = 0 ==========

tp = Array{Vector{Float64}}(nrow(df_portfolio))
tp_0 = zeros(Float64, nrow(df_portfolio))

for d = 1:nrow(df_portfolio)
  # d = 4
  tp[d] = zeros(Float64, T)
  prob_sx =
    convert(Array, prob[:sx]) * df_portfolio[d, :sx_be_fac]
  τ = t_0 - df_portfolio[d, :t_start]
  for t = (T-1-τ):-1:(1)
    tp[d][t] = tpberec(tp[d][t + 1], t, τ, prob_sx)
  end
  tp_0[d] = tpberec(tp[d][1], 0, τ, prob_sx)
end

for d = 1:nrow(df_portfolio)
  for t = 1:(T + df_portfolio[d, :t_start])
    @test VORMI.tpg(t, cap_mkt.rfr.x, liab_ins.mps[d]) ≈
          df_portfolio[d, :n] * tp[d][t]
  end
end

lx = zeros(Float64, T+1)
tp_all_0 = 0.0
for d = 1:nrow(df_portfolio)
  τ = t_0 - df_portfolio[d, :t_start]
  lx[1] = 1.0
  tp_all_0 += lx[1] * tp_0[d] * df_portfolio[d, :n]
  for t = 1:T
    if t + τ - t_0 <= T
      lx[t + 1] = lx[t] * prob[t + τ - t_0, :px]
    end
  end
end

@test tp_all_0 ≈ proj.val_0[1, :tpg]

## Cashflows year 1 =============================================

## bi-quotient --------------------------------------------------
@test VORMI.bonusrate(1, y_invest[1], liab_ins.mps[1], dyn) ≈
      t_bonus_quota[1]
@test VORMI.yield(1, cap_mkt.rfr) ≈ rfr[1]
@test VORMI.yield(1, cap_mkt.stock) ≈ y_stock[1]
ind_bonus_1 =
  y_stock[1] / (t_bonus_quota[1] + liab_ins.mps[1].rfr_price[1])

@test VORMI.yield(0, cap_mkt.stock) ≈ cap_mkt.stock.yield_0
@test VORMI.yield(0, cap_mkt.rfr) ≈ cap_mkt.rfr.yield_0
ind_bonus_hypo =
  cap_mkt.stock.yield_0 /
  (liab_ins.mps[1].bonus_rate_hypo + liab_ins.mps[1].rfr_price_0)

bi_quot_1 =ind_bonus_1 / ind_bonus_hypo

@test bi_quot_1 ≈ VORMI.biquotient(1, y_invest[1], cap_mkt,
                                   invs, liab_ins.mps[1], dyn)

## In the text we assume that b^C,hypo does not depend on C
for d = 2:nrow(df_portfolio)
  @test df_portfolio[1,:bonus_rate_hypo] ≈
        df_portfolio[d,:bonus_rate_hypo]
end

## δ_sx_one, qx_one, sx_one, px_one are vectors
## over all model points for time t == 1
δ_sx_one =
  Float64[VORMI.δsx(1, cap_mkt, invs, liab_ins.mps[d], dyn)
   for d = 1:nrow(df_portfolio)]

t = 1
for d = 1:nrow(df_portfolio)
  if length(sx_basis[d]) >= t
    @test sx_basis[d][t] ≈ liab_ins.mps[d].prob[t, :sx]
  end
end

sx_one = δ_sx_one .* Float64[sx_basis[d][1] for d = 1:nrow(df_portfolio)]

qx_one =
  Float64[liab_ins.mps[d].prob[t, :qx]
          for d = 1:nrow(df_portfolio)]
px_one = 1 .- qx_one .- sx_one

## Cashflows ----------------------------------------------------

cf_prem_one =
  sum([liab_ins.mps[d].β[1,:prem] for d = 1:nrow(df_portfolio)])
cf_λ_boy_one =
  -sum([liab_ins.mps[d].λ[1,:boy] for d = 1:nrow(df_portfolio)])
cf_λ_eoy_one =
  -sum([liab_ins.mps[d].λ[1,:eoy] *
          (1 + liab_ins.mps[d].λ[1,:infl])
        for d = 1:nrow(df_portfolio)]) -
  invs.igs[:IGCash].mv_0 *
  (1 + (cf_prem_one + cf_λ_boy_one) / invs.mv_0) *
  invs.igs[:IGCash].cost.rel[1] *
  (1 + invs.igs[:IGCash].cost.infl_rel[1]) -
  invs.igs[:IGStock].mv_0 *
  (1 + (cf_prem_one + cf_λ_boy_one) / invs.mv_0) *
  invs.igs[:IGStock].cost.rel[1] *
  (1 + invs.igs[:IGStock].cost.infl_rel[1]) -
  invs.igs[:IGCash].cost.abs[1] *
  (1 + invs.igs[:IGCash].cost.infl_abs[1]) -
  invs.igs[:IGStock].cost.abs[1] *
  (1 + invs.igs[:IGStock].cost.infl_abs[1])
cf_qx_one = -sum([qx_one[d] *liab_ins.mps[d].β[1,:qx]
                  for d = 1:nrow(df_portfolio)])
cf_sx_one = -sum([sx_one[d] *liab_ins.mps[d].β[1,:sx]
                  for d = 1:nrow(df_portfolio)])
cf_px_one = -sum([px_one[d] *liab_ins.mps[d].β[1,:px]
                  for d = 1:nrow(df_portfolio)])
cf_bonus_one =
  -sum([t_bonus_quota[1] * liab_ins.mps[d].tpg_price_0
        for d = 1:5])
cf_invest_one =
  (invs.mv_0+cf_prem_one + cf_λ_boy_one) * y_invest[1]

@test cf_prem_one ≈ ins_sum * sum(df_portfolio[:,:n]) *
                    prem_price_ratio
@test cf_qx_one ≈ -sum([qx_one[d] * df_portfolio[d,:n] * ins_sum
                        for d = 1:nrow(df_portfolio)])
@test cf_sx_one ≈ -sum([sx_one[d] * sx_fac * df_portfolio[d,:n] *
                        (T - d +1) *  ins_sum * prem_price_ratio
                        for d = 1:nrow(df_portfolio)])
@test cf_px_one ≈ -px_one[1] * df_portfolio[1,:n] *  ins_sum
@test cf_bonus_one ≈ proj.cf[1, :bonus]
@test cf_prem_one ≈ proj.cf[t, :prem]
@test cf_λ_boy_one ≈ proj.cf[t, :λ_boy]
@test cf_λ_eoy_one ≈ proj.cf[t, :λ_eoy]
@test cf_qx_one ≈ proj.cf[t, :qx]
@test cf_sx_one ≈ proj.cf[t, :sx]
@test cf_px_one ≈ proj.cf[t, :px]
@test cf_invest_one ≈ proj.cf[t, :invest]
@test cf_l_other_total[1] ≈ proj.cf[t, :l_other]


## technical provisions for guaranteed benefits -----------------

## the following is used in the text
@test(abs(δ_sx_one[1]-1) > eps(1.) ? true : false)

probs = Array{DataFrame}(nrow(df_portfolio))
for d = 1:nrow(df_portfolio)
  probs[d] = deepcopy(prob_be)
  probs[d][:sx] =
    δ_sx_one[d] * convert(Array, probs[d][:sx]) *
    df_portfolio[d, :sx_be_fac]
  probs[d][:px] = 1 .- probs[d][:qx] - probs[d][:sx]
end

## We recalculate technical provisions with updated sx
for d = 1:nrow(df_portfolio)
  tp[d] = zeros(Float64, T)
  τ = t_0 - df_portfolio[d, :t_start]
  for t = (T-1-τ):-1:(1)
    tp[d][t] = tpberec(tp[d][t + 1], t, τ, probs[d][:sx])
  end
end

tpg_0_1_1 =
  sum([probs[d][T-d+1, :px] * tp[d][1] *
         df_portfolio[d, :n]
       for d = 1:nrow(df_portfolio)])

@test tpg_0_1_1 ≈ proj.val[1,:tpg]
@test proj.val[1,:tpg]-proj.val_0[1,:tpg] ≈ -proj.cf[1,:Δtpg]

# Profit and tax ------------------------------------------------
@test sum(convert(Array, proj.cf[1, [:prem, :λ_boy, :λ_eoy, :qx,
                                     :sx, :px, :invest, :Δtpg, :l_other, :bonus]])) ≈
      proj.cf[1, :profit]

tax_pre = zeros(Float64, T)
tax= zeros(Float64, T)
tax_credit = tax_credit_0
tax_credit_vec = zeros(Float64, T)

for t = 1:T
  tax_pre[t] = tax_rate * proj.cf[t,:profit]
  if tax_credit > tax_pre[t]
    tax[t] = 0
    tax_credit -= tax_pre[t]
  else
    tax[t] = tax_pre[t] - tax_credit
    tax_credit = 0
  end
  tax_credit_vec[t] = tax_credit
end

@test tax ≈ -proj.cf[:,:tax]

## Dividends ----------------------------------------------------
surp_quota = zeros(Float64,T)
for t = 1:T
  surp_quota[t] =
    proj.val[t, :invest] /
    (proj.val[t, :tpg] + proj.val[t, :l_other]) -1
end

invest_eoy_prev =
  Float64[t == 1 ?
            proj.val_0[t, :invest] :
            proj.val[t-1, :invest] for t = 1:T]
invest_boy =
  convert(Array,
          invest_eoy_prev + proj.cf[:prem] + proj.cf[:λ_boy])
invest_eoy_pre_divid =
  invest_eoy_prev + proj.cf[:profit] + proj.cf[:tax] -
  proj.cf[:Δtpg] + proj.cf[:gc]

@test invest_boy ≈ invs.mv_boy
@test invest_eoy_pre_divid ≈
      Float64[VORMI.investpredivid(t, invs, proj) for t = 1:T]

val_liab = convert(Array, proj.val[:tpg] .- proj.val[:l_other])
q_surp =
  ((invest_eoy_pre_divid .- proj.val[:tpg] .-
    proj.val[:l_other]) ./
  (proj.val[:tpg] + proj.val[:l_other]))

## The following was used in the text:
@test q_surp[1] > dyn.quota_surp
@test min.(0, convert(Array, (1+dyn.quota_surp) *
      (proj.val[:tpg] + proj.val[:l_other]) - invest_eoy_pre_divid)) ≈
      proj.cf[:divid]

## Dividend mechanism works:
for t = 1:T
  if (proj.cf[t, :divid] < 0) &
      (proj.val[t, :tpg] + proj.val[t, :l_other] > 0)
    @test dyn.quota_surp ≈
          (proj.val[t, :invest] -
          proj.val[t, :tpg] - proj.val[t, :l_other]) / (proj.val[t, :tpg] + proj.val[t, :l_other])
  end
  if proj.cf[t, :divid] ≥ 0
    @test proj.cf[t, :divid] ≈ 0
  end
end

@test proj.val[1, :invest] ≈
      invest_eoy_pre_divid[1] + proj.cf[1, :divid]

## Balance sheet for ≥ 0 ========================================

fdb = zeros(Float64, T)
for t = (T-1):-1:1
  fdb[t] =
    (fdb[t + 1] -  proj.cf[t + 1, :bonus]) /
    (1 .+ rfr[t + 1])
end
fdb_0 =
  (fdb[1] -  proj.cf[1, :bonus]) /  (1 .+ rfr[1])

@test -cumprod(1 ./ (1 .+ rfr)) ⋅ proj.cf[:bonus] ≈
      proj.val_0[1,:bonus]
@test fdb_0 ≈ proj.val_0[1,:bonus]
@test fdb ≈ proj.val[:bonus]

balance = vcat(proj.val_0, proj.val)
## here reserves for boni are considered part of capital.
@test balance[:surplus] ≈
      balance[:invest] - balance[:tpg] - balance[:l_other]

## recall that balance[t, :] = proj.val[t-1, :] for t>1
for t in 2:(T+1)
  for w in [:invest, :tpg, :l_other, :surplus, :bonus, :cost_prov]
     @test balance[t, w] ≈ proj.val[t-1, w]
  end
end

## technical provisions for absolute costs & investment costs ---

# provisions absolute costs & investment costs are correct:
for τ = 1:5
  x = balance[τ, :cost_prov]
  for t = τ:T
    x *= (1 + rfr[t])
    x -= sum(convert(Array,
                     balance[t,
                             [:tpg, :bonus,
                              :l_other, :cost_prov]])) *
      invs.igs[:IGCash].cost.cum_infl_rel[t] *
      invs.igs[:IGCash].cost.rel[t]
    x -= proj.fixed_cost_gc[t]
  end
  @test x ≈ 0 atol=1.0e-14
end


tpgprev(tp_next, t) =
  (tp_next + cost_abs[t + 1])/(1 + rfr[t + 1])

tp_cost_abs = zeros(Float64, T)
for t = (T-1):-1:1
  tp_cost_abs[t] = tpgprev(tp_cost_abs[t + 1], t)
end
tp_cost_abs_0 = tpgprev(tp_cost_abs[1], 0)

@test [tp_cost_abs_0; tp_cost_abs] ≈
      Float64[VORMI.tpgfixed(t, cap_mkt.rfr.x[1:liab_ins.dur],
              proj.fixed_cost_gc) for t in 0:T]

## S2 Example ###################################################

liabs_mod_0 = balance[1,:tpg] +balance[1,:bonus] + balance[1, :cost_prov]
assets_mod_0 = balance[1,:invest] + proj.tax_credit_0
bof_0 = assets_mod_0 - liabs_mod_0
symb_bal = [:invest, :tpg, :l_other, :surplus, :bonus]
@test convert(Array, balance[1,symb_bal]) ≈
      convert(Array, s2.balance[1, symb_bal])
@test VORMI.bof(s2, :be) ≈ bof_0

## S2 Example Interest ------------------------------------------

ind_mkt = findin(ds2[:mdl], [:S2Mkt])[1]
ind_mkt_int = findin(ds2_mkt[:mdl], [:S2MktInt])[1]
s2_mkt_int = s2.mds[ind_mkt].mds[ind_mkt_int]
rfr_up = VORMI.rfrshock(cap_mkt.rfr.x,
                        s2_mkt_int,
                        :spot_up)
rfr_down = VORMI.rfrshock(cap_mkt.rfr.x,
                          s2_mkt_int,
                          :spot_down)

@test rfr ≈ VORMI.spot2forw(VORMI.forw2spot(rfr))
@test rfr_down ≈
      VORMI.spot2forw(VORMI.forw2spot(rfr) .*
      (1 .+ ds2_mkt_int[:shock][:spot_down][1:T]))
@test rfr_up ≈ VORMI.spot2forw(
                  VORMI.forw2spot(rfr) .+
                  max.(ds2_mkt_int[:spot_up_abs_min],
                      VORMI.forw2spot(rfr) .*
                      ds2_mkt_int[:shock][:spot_up][1:T]))

## S2 Example Equity --------------------------------------------
ind_mkt = findin( ds2[:mdl], [:S2Mkt])[1]
ind_mkt_eq = findin(ds2_mkt[:mdl], [:S2MktEq])[1]
s2_mkt_eq = s2.mds[ind_mkt].mds[ind_mkt_eq]
bal = s2_mkt_eq.balance
@test bal[bal[:scen] .== :type_1,:invest][1] ≈
      (1 + eq_shock[:type_1]) * sum(df_stock[:mv_0]) + sum(df_cash[:mv_0])

## S2 Example Market risk----------------------------------------
ind_mkt = findin( ds2[:mdl], [:S2Mkt])[1]
s2_mkt = s2.mds[ind_mkt]
@test s2_mkt_int.scen_up == false
corr_mkt = s2_mkt.corr_down[1:2,1:2]
scr_mkt_net = [s2_mkt_int.scr[NET], s2_mkt_eq.scr[NET]]
scr_mkt_gross = [s2_mkt_int.scr[GROSS], s2_mkt_eq.scr[GROSS]]
@test sqrt(scr_mkt_net ⋅ (corr_mkt * scr_mkt_net)) ≈
      s2_mkt.scr[NET]
@test sqrt(scr_mkt_gross ⋅ (corr_mkt * scr_mkt_gross)) ≈
      s2_mkt.scr[GROSS]

## Default Risk type 1 -----------------------------------------
ind_def = findin( ds2[:mdl], [:S2Def])[1]
s2_def = s2.mds[ind_def]
inv_len = length(invs.igs[:IGCash].investments)
@test inv_len == 2

accs_mv_0 =
  Float64[invs.igs[:IGCash].investments[i].mv_0
          for i = 1:inv_len]
accs_cqs =
  Int[parse(Int, string(string(invs.igs[:IGCash].
                                 investments[i].cqs)[end]))
      for i in 1:inv_len]
accs_tlgd =
  Float64[s2_def.mds[1].tlgd[accs_cqs[i]+1] for i in 1:inv_len]
accs_slgd =
  Float64[s2_def.mds[1].slgd[accs_cqs[i]+1] for i in 1:inv_len]
accs_defu =
  Float64[s2_def.mds[1].u[accs_cqs[i]+1,accs_cqs[j]+1]
          for i in 1:inv_len, j in 1:inv_len]
accs_defv =
  Float64[s2_def.mds[1].v[accs_cqs[i]+1]
          for i in 1:inv_len]

accs_var_t = accs_tlgd ⋅ (accs_defu * accs_tlgd)
accs_var_s = accs_slgd ⋅ accs_defv
accs_var = accs_var_t + accs_var_s
accs_sigma_norm = sqrt(accs_var) / sum(accs_tlgd)
accs_scr_low = s2.mds[ind_def].mds[1].scr_par[:low][1]
accs_scr_low_fac = s2.mds[ind_def].mds[1].scr_par[:low][2]
accs_scr_def = accs_scr_low_fac * sqrt(accs_var)

@test s2_def.scr[NET] ≈ accs_scr_def
@test s2_def.scr[GROSS] ≈ accs_scr_def

## Life mortality -----------------------------------------------
ind_life = findin( ds2[:mdl], [:S2Life])[1]
ind_life_qx = 1
s2_life_qx = s2.mds[ind_life].mds[ind_life_qx]

@test collect(keys(s2_life_qx.shock))[1] == :qx
@test s2_life_qx.shock[:qx] > 0
@test ! s2_life_qx.mp_select[:qx][1]
for i = 2:length(s2_life_qx.mp_select[:qx])
  @test s2_life_qx.mp_select[:qx][i]
end

## Life longevity -----------------------------------------------

ind_life_px = 2
s2_life_px = s2.mds[ind_life].mds[ind_life_px]
@test collect(keys(s2_life_px.shock))[1] == :px
@test s2_life_px.shock[:px] < 0
for i = 2:length(s2_life_px.mp_select[:px])
  @test !s2_life_px.mp_select[:px][i]
end

## Life surrender -----------------------------------------------

ind_life_sx = 4
s2_life_sx = s2.mds[ind_life].mds[ind_life_sx]

@test (sort(collect(keys(s2_life_sx.shock))) ==
         sort([:sx_down, :sx_up,
               :sx_mass_other, :sx_mass_pension]))

for m in 2: length(s2_life_sx.mp_select[:sx_down])
  @test ! s2_life_sx.mp_select[:sx_down][m]
end
for m in 2:length(s2_life_sx.mp_select[:sx_up])
  @test s2_life_sx.mp_select[:sx_up][m]
end
for m in 2:length(s2_life_sx.mp_select[:sx_mass_other])
  @test s2_life_sx.mp_select[:sx_mass_other][m]
end

## Life cost ----------------------------------------------------

ind_life_cost = 5
s2_life_cost = s2.mds[ind_life].mds[ind_life_cost]
@test collect(keys(s2_life_cost.shock)) == [:cost]

## Life cat -----------------------------------------------------
s2.mds[ind_life].mds
ind_life_cat = 7
s2_life_cat = s2.mds[ind_life].mds[ind_life_cat]
@test collect(keys(s2_life_cat.shock)) == [:cat]
for m = 2:length(s2_life_cat.mp_select[:cat])
  @test s2_life_cat.mp_select[:cat][m]
end

## Life aggregation ---------------------------------------------

s2_life = s2.mds[ind_life]

ind_life_risks =
  [ind_life_qx, ind_life_sx, ind_life_cost, ind_life_cat]
life_corr = s2_life.corr[ind_life_risks, ind_life_risks]

scrs_life_net = [s2_life_qx.scr[NET],
                 s2_life_sx.scr[NET],
                 s2_life_cost.scr[NET],
                 s2_life_cat.scr[NET]]
scrs_life = [s2_life_qx.scr[GROSS],
             s2_life_sx.scr[GROSS],
             s2_life_cost.scr[GROSS],
             s2_life_cat.scr[GROSS]]

@test sqrt(scrs_life ⋅ (life_corr * scrs_life)) ≈
      s2_life.scr[GROSS]
@test sqrt(scrs_life_net ⋅ (life_corr * scrs_life_net)) ≈
      s2_life.scr[NET]

## BSCR =========================================================

ind = [ind_mkt, ind_def, ind_life]

bscr_corr = ds2[:corr][ind, ind]
bscrs_gross =
  [s2_mkt.scr[GROSS], s2_def.scr[GROSS], s2_life.scr[GROSS]]

bscrs_net = [s2_mkt.scr[NET], s2_def.scr[NET], s2_life.scr[NET]]

@test sqrt(bscrs_net ⋅ (bscr_corr * bscrs_net)) ≈ s2.bscr[NET]
@test sqrt(bscrs_gross ⋅ (bscr_corr * bscrs_gross)) ≈
      s2.bscr[GROSS]

## operational risk =============================================
@test s2.op.comp_tp ≈
      (proj.val_0[1, :tpg] + proj.val_0[1, :bonus]) * s2.op.fac[:tp]
@test s2.op.comp_prem  ≈
      s2.op.fac[:prem] * s2.op.prem_earned + s2.op.fac[:prem] *  max(0, s2.op.prem_earned -
                s2.op.fac[:prem_py] * s2.op.prem_earned_prev)

## Adjustments ==================================================
## Adj technical provisions -------------------------------------
@test -max(0.0,  min(s2.bscr[GROSS] - s2.bscr[NET],
                     VORMI.fdb(s2, :be))) ≈
      s2.adj_tp
@test s2.adj_dt == 0
@test s2.liabs_mod ≈
      s2.balance[1,:tpg] + s2.balance[1,:bonus] + s2.balance[1,:cost_prov]
@test s2.invest_mod ≈ s2.balance[1,:invest]

ds2[:coc]

coc = 0.06
balance = vcat(proj.val_0, proj.val)
length(balance)
tpbe =
  convert(Array,
          balance[:tpg] + balance[:bonus] + balance[:cost_prov])
src_future = (tpbe * s2.scr / tpbe[1])[1:T]
discount = 1 ./ cumprod(1 .+ rfr)
risk_margin = coc * src_future ⋅ discount

@test s2.risk_margin ≈ risk_margin

@test s2.scr_ratio ≈
      (s2.invest_mod - s2.liabs_mod - s2.risk_margin) / s2.scr
@test s2.scr_ratio > 1

#################################################################
println("End S2Life test")

sqrt
