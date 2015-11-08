## S2MktInt -----------------------------------------------------
function S2MktInt(ds2_mkt_int::Dict{Symbol, Any})
  shock_object = :CapMkt
  shock = ds2_mkt_int[:shock]
  spot_up_abs_min = ds2_mkt_int[:spot_up_abs_min]
  balance = DataFrame()
  scr = zeros(Float64, 2)
  scen_up = false
  return S2MktInt(shock_object,
                  shock,
                  spot_up_abs_min, balance, scr, scen_up)
end

function S2MktInt(param::ProjParam,
                  s2_balance::DataFrame,
                  ds2_mkt_int::Dict{Symbol, Any})
  p = deepcopy(param)
  mkt_int = S2MktInt(ds2_mkt_int)
  mkt_int.balance = deepcopy(s2_balance)
  for int_type_symb in collect(keys(mkt_int.shock))
    append!(mkt_int.balance,
            s2bal(p, mkt_int,
                  (inv, s2_int) ->
                  mktintshock!(inv, s2_int, int_type_symb),
                  int_type_symb))
  end
  scr!(mkt_int)
  return mkt_int
end

## S2MktEq ------------------------------------------------------
function S2MktEq(ds2_mkt_eq::Dict{Symbol, Any})
  shock_object = :CapMkt_AdjVal0
  eq_type = Dict{Symbol, Symbol}()
  balance = DataFrame()
  corr = ds2_mkt_eq[:corr]
  shock = ds2_mkt_eq[:shock]
  scr = zeros(Float64, 2)
  return S2MktEq(shock_object, eq_type, shock,
                 balance, corr, scr)
end

function S2MktEq(param::ProjParam,
                 s2_balance::DataFrame,
                 ds2_mkt_eq::Dict{Symbol, Any},
                 eq2type)
  p = deepcopy(param)
  mkt_eq = S2MktEq(ds2_mkt_eq)
  merge!(mkt_eq.eq2type, eq2type)
  mkt_eq.balance = deepcopy(s2_balance)
  for shock_symb in collect(keys(mkt_eq.shock))
    append!(mkt_eq.balance,
            s2bal(p, mkt_eq,
                  (invs, s2_eq) ->
                  mkteqshock!(invs, s2_eq, shock_symb),
                  shock_symb))
  end
  scr!(mkt_eq)
  return mkt_eq
end

## S2Mkt --------------------------------------------------------
function S2Mkt(ds2_mkt::Dict{Symbol, Any})
  mds = Array(S2Module, 0)
  corr_up = ds2_mkt[:corr](ds2_mkt[:raw], ds2_mkt[:adj], :up)
  corr_down = ds2_mkt[:corr](ds2_mkt[:raw], ds2_mkt[:adj], :down)
  scr = zeros(Float64, 2)
  return S2Mkt(mds, corr_up, corr_down, scr)
end


function S2Mkt(param::ProjParam,
               s2_balance::DataFrame,
               eq2type::Dict,
               ds2_mkt_all::Dict)
  p = deepcopy(param)
  mkt = S2Mkt(ds2_mkt_all[:mkt])
  push!(mkt.mds, S2MktInt(p, s2_balance,
                 ds2_mkt_all[:mkt_int]))
  push!(mkt.mds, S2MktEq(p, s2_balance,
                 ds2_mkt_all[:mkt_eq], eq2type))
  push!(mkt.mds, S2MktProp(p, s2_balance))
  push!(mkt.mds, S2MktSpread(p, s2_balance))
  push!(mkt.mds, S2MktFx(p, s2_balance))
  push!(mkt.mds, S2MktConc(p, s2_balance))
  scen_up = false
  for i = 1:length(mkt.mds)
    if :scen_up in fieldnames(mkt.mds[i])
      scen_up = mkt.mds[i].scen_up
    end
  end
  corr = (scen_up ? mkt.corr_up : mkt.corr_down)
  mkt.scr = scr(mkt, corr)
  return mkt
end

## S2Def1 -------------------------------------------------------
function S2Def1(ds2_def_1)
  tlgd = Array(Float64, 0)
  slgd = Array(Float64, 0)
  u = Array(Float64, 0, 0)
  v = Array(Float64, 0)
  scr_par = Dict{Symbol, Vector{Float64}}()
  for i = 1:nrow(ds2_def_1[:scr_par])
    merge!(scr_par,
           Dict(ds2_def_1[:scr_par][i, :range] =>
            [ds2_def_1[:scr_par][i, :threshold_upper],
             ds2_def_1[:scr_par][i, :multiplier]]))
  end
  scr = zeros(Float64, 2)
  return S2Def1(tlgd, slgd, u, v, scr_par, scr)
end

function S2Def1(param::ProjParam,
                ds2_def_1::Dict{Symbol, Any})
  p = deepcopy(param)
  def = S2Def1(ds2_def_1)
  cqs_vec = filter(x -> ismatch(r"cqs", string(x)),
                   names(ds2_def_1[:prob]))
  prob = [ds2_def_1[:prob][1, cqs] for cqs in cqs_vec]
  def.tlgd = zeros(Float64, length(cqs_vec))
  def.slgd = zeros(Float64, length(cqs_vec))
  def.u = Array(Float64, length(cqs_vec), length(cqs_vec))
  def.v = Array(Float64, length(cqs_vec))

  def.v = 1.5 * prob .* (1 .- prob) ./ (2.5 .- prob)
  for i = 1:size(def.u,1), j = 1:1:size(def.u,2)
    def.u[i,j] =
      (1-prob[i]) * prob[i] * (1-prob[j]) * prob[j] /
      (1.25 * (prob[i] + prob[j]) - prob[i] * prob[j])
  end
  invs = InvPort(p.t_0, p.T, p.cap_mkt, p.invs_par...)
  for i = 1:length(invs.igs[:IGCash].investments)
    j = indexin([invs.igs[:IGCash].investments[i].cqs],
                cqs_vec)[1]
    lgd =
      invs.igs[:IGCash].investments[i].lgd *
      invs.igs[:IGCash].investments[i].mv_0
    def.tlgd[j] += lgd
    def.slgd[j] += lgd * lgd
  end
  scr!(def)
  return def
end

## S2Def --------------------------------------------------------
function S2Def(ds2_def)
  mds = Array(S2Module, 0)
  corr = ds2_def[:corr]
  scr = zeros(Float64, 2)
  return S2Def(mds, corr, scr)
end

function S2Def(param::ProjParam,
               s2_balance::DataFrame,
               ds2_def_all::Dict)
  p = deepcopy(param)
  def = S2Def(ds2_def_all[:def])
  push!(def.mds, S2Def1(p, ds2_def_all[:def_1]))
  push!(def.mds, S2Def2(p, s2_balance))
  def.scr = scr(def, def.corr)
  return def
end

## S2LifeBio (mortality risk, longevity risk, lapse risk---------
function S2LifeBio(ds2_bio::Dict{Symbol, Any})
  shock_object = :LiabIns
  if :shock_param in collect(keys(ds2_bio))
    shock_param = ds2_bio[:shock_param]
  else
    shock_param = Dict()
  end
  shock = ds2_bio[:shock]
  balance = DataFrame()
  mp_select = Dict{Symbol, Vector{Bool}}()
  scr = zeros(Float64, 2)
  return S2LifeBio(shock_object, shock,
                   shock_param,
                   balance, mp_select, scr)
end

function S2LifeBio(param::ProjParam,
                   s2_balance::DataFrame,
                   ds2_bio::Dict{Symbol, Any})
  p = deepcopy(param)
  bio = S2LifeBio(ds2_bio)
  bio.balance = deepcopy(s2_balance)
  select!(p, bio)
  for symb in collect(keys(bio.shock))
    append!(bio.balance,
            s2bal(p, bio,
                  (l_ins, wx) -> bioshock!(l_ins, wx, symb),
                  symb))
  end
  scr!(bio)
  return bio
end

## S2LifeCost ---------------------------------------------------
function S2LifeCost(ds2_cost::Dict{Symbol, Any})
  shock_object = :InvPort_LiabIns
  shock = ds2_cost[:shock]
  if :shock_param in collect(keys(ds2_cost))
    shock_param = ds2_cost[:shock_param]
  else
    shock_param = Dict()
  end
  balance = DataFrame()
  scr = zeros(Float64, 2)
  return S2LifeCost(shock_object, shock, shock_param,
                    balance, scr)
end

function S2LifeCost(param::ProjParam,
                    s2_balance::DataFrame,
                    ds2_cost::Dict{Symbol, Any})
  p = deepcopy(param)
  cost = S2LifeCost(ds2_cost)
  cost.balance = deepcopy(s2_balance)
  for symb in collect(keys(cost.shock))
    append!(cost.balance,
            s2bal(p, cost, (invs, l_ins, cst) ->
                  costshock!(invs, l_ins, cst), symb))
  end
  scr!(cost)
  return cost
end

## S2Life -------------------------------------------------------
function S2Life(ds2_life::Dict{Symbol, Any})
  mds = Array(S2Module, 0)
  corr = ds2_life[:corr]
  scr = zeros(Float64, 2)
  return S2Life(mds, corr, scr)
end

function S2Life(param::ProjParam,
                s2_balance::DataFrame,
                d_life::Dict)
  p = deepcopy(param)
  life = S2Life(d_life[:life])
  push!(life.mds, S2LifeBio(p, s2_balance, d_life[:life_qx]))
  push!(life.mds, S2LifeBio(p, s2_balance, d_life[:life_px]))
  push!(life.mds, S2LifeMorb(p, s2_balance))
  push!(life.mds, S2LifeBio(p, s2_balance, d_life[:life_sx]))
  push!(life.mds, S2LifeCost(p, s2_balance, d_life[:life_cost]))
  push!(life.mds, S2LifeRevision(p, s2_balance))
  push!(life.mds, S2LifeBio(p, s2_balance, d_life[:life_cat]))
  life.scr = scr(life, life.corr)
  return life
end

## S2Op ---------------------------------------------------------
S2Op(ds2_op::Dict{Symbol, Float64}, s2_op::Dict) =
  S2Op(ds2_op,
       s2_op[:prem_earned],
       s2_op[:prem_earned_prev],
       s2_op[:tp], 0, 0, s2_op[:cost_ul], 0.0)

## S2 -----------------------------------------------------------
function   S2(ds2_op, s2_op)
  mds = Array(S2Module, 0)
  balance = DataFrame()
  corr = zeros(Float64, 5, 5)
  bscr = zeros(Float64, 2)
  adj_tp = 0.0
  adj_dt = 0.0
  op = S2Op(ds2_op, s2_op)
  scr = 0.0
  S2(mds, balance, corr, bscr, adj_tp, adj_dt, op, scr, 0, 0, 0)
end

function  S2(param::ProjParam,
             eq2type::Dict{Symbol, Symbol},
             ds2_mkt_all::Dict,
             ds2_def_all::Dict,
             ds2_life_all::Dict,
             s2_op::Dict,
             ds2_op::Dict{Symbol, Float64},
             ds2::Dict{Symbol, Any})
  p = deepcopy(param)
  s2 = S2(ds2_op, s2_op)
  s2.corr = ds2[:corr]
  s2.balance = s2bal(p)
  s2.op.tp = s2.balance[1, :tpg] + s2.balance[1, :bonus]
  push!(s2.mds, S2Mkt(p, s2.balance, eq2type, ds2_mkt_all))
  push!(s2.mds, S2Def(p, s2.balance, ds2_def_all))
  push!(s2.mds, S2Life(p, s2.balance, ds2_life_all))
  push!(s2.mds, S2Health(p, s2.balance))
  push!(s2.mds, S2NonLife(p, s2.balance))
  scr!(s2, p.tax_credit_0)
  return s2
end
