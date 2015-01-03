
## constructors =================================================
## S2MktInt -----------------------------------------------------
function S2MktInt(ds2_mkt_int::Dict{Symbol, Any})
  shock_object = :CapMkt
  shock_type = collect(keys(ds2_mkt_int[:shock]))
  shock = ds2_mkt_int[:shock]
  spot_up_abs_min = ds2_mkt_int[:spot_up_abs_min]
  balance = DataFrame()
  scr = zeros(Float64, 2)
  scen_up = false
  return S2MktInt(shock_object, shock_type, shock,
                  spot_up_abs_min, balance, scr, scen_up)
end

function S2MktInt(p::ProjParam,
                  s2_balance::DataFrame,
                  ds2_mkt_int::Dict{Symbol, Any})
  s2_mkt_int = S2MktInt(ds2_mkt_int)
  s2_mkt_int.balance = deepcopy(s2_balance)
  for int_type_symb in s2_mkt_int.shock_type
  end
  for int_type_symb in s2_mkt_int.shock_type
    append!(s2_mkt_int.balance,
            s2bal(p, s2_mkt_int,
                  (inv, s2_int) ->
                  mktintshock!(inv, s2_int, int_type_symb),
                  int_type_symb))
  end
  scr!(s2_mkt_int)
  return s2_mkt_int
end

## S2MktEq ------------------------------------------------------
function S2MktEq(ds2_mkt_eq::Dict{Symbol, Any})
  shock_object = :InvPort
  shock_type = collect(keys(ds2_mkt_eq[:shock]))
  eq_type = Dict{Symbol, Symbol}()
  balance = DataFrame()
  corr = ds2_mkt_eq[:corr]
  shock = ds2_mkt_eq[:shock]
  scr = zeros(Float64, 2)
  return S2MktEq(shock_object, shock_type, eq_type, shock,
                 balance, corr, scr)
end

function S2MktEq(p::ProjParam,
                 s2_balance::DataFrame,
                 s2_mkt_eq::Dict{Symbol, Any},
                 eq2type)
  s2mkt_eq = S2MktEq(s2_mkt_eq)
  merge!(s2mkt_eq.eq2type, eq2type)
  s2mkt_eq.balance = deepcopy(s2_balance)
  for eq_type_symb in s2mkt_eq.shock_type
    append!(s2mkt_eq.balance,
            s2bal(p, s2mkt_eq,
                  (inv, s2_eq) ->
                  mkteqshock!(inv, s2_eq, eq_type_symb),
                  eq_type_symb))
  end
  scr!(s2mkt_eq)
  return s2mkt_eq
end


## S2Mkt --------------------------------------------------------
function S2Mkt(s2_mkt::Dict{Symbol, Any})
  mds_symb = [:S2MktInt, :S2MktEq, :S2MktSpread]
  mds = Array(S2Module, 0)
  mds_scr = zeros(Float64, length(mds_symb), 2)
  corr_up = s2_mkt[:corr](s2_mkt[:raw], s2_mkt[:adj], :up)
  corr_down = s2_mkt[:corr](s2_mkt[:raw], s2_mkt[:adj], :down)
  scr = zeros(Float64, 2)
  return S2Mkt(mds_symb, mds, mds_scr, corr_up, corr_down, scr)
end

function S2Mkt(p::ProjParam,
               s2_balance::DataFrame,
               ds2_mkt_int::Dict{Symbol, Any},
               ds2_mkt_eq::Dict{Symbol, Any},
               eq2type::Dict{Symbol, Symbol},
               ds2_mkt::Dict{Symbol, Any} )
  s2mkt = S2Mkt(ds2_mkt)
  push!(s2mkt.mds, S2MktInt(p, s2_balance, ds2_mkt_int))
  push!(s2mkt.mds, S2MktEq(p, s2_balance, ds2_mkt_eq, eq2type))
  push!(s2mkt.mds, S2MktProp(p, s2_balance))
  push!(s2mkt.mds, S2MktSpread(p, s2_balance))
  push!(s2mkt.mds, S2MktFx(p, s2_balance))
  push!(s2mkt.mds, S2MktConc(p, s2_balance))
  scr!(s2mkt)
  return s2mkt
end

## S2Def1 -------------------------------------------------------
function S2Def1(ds2_def_1)
  tlgd = Array(Float64, 0)
  slgd = Array(Float64, 0)
  u = Array(Float64, 0, 0)
  v = Array(Float64, 0)
  scr = zeros(Float64, 2)
  return S2Def1(tlgd, slgd, u, v, scr)
end

function S2Def1(p::ProjParam,
                ds2_def_1::Dict{Symbol, Any})
  cqs = filter(x -> ismatch(r"cqs", string(x)),
               names(ds2_def_1[:prob]))
end

## S2Def --------------------------------------------------------
function S2Def(p::ProjParam,
               s2_balance::DataFrame,
               ds2_def_1::Dict{Symbol, Any},
               ds2_def::Dict{Symbol, Any})
  s2_def = S2Def(ds2_def)
  push!(s2_def1.mds, S2Def1(p, ds2_def_1))
  push!(s2_def2.mds, S2Def2(p, s2_balance))
  scr!(s2_def)
  return s2_def
end

## S2Op ---------------------------------------------------------
S2Op(s2_op::Dict{Symbol, Float64}) =
  S2Op(s2_op[:prem_earned], s2_op[:prem_earned_prev],
       s2_op[:tp], s2_op[:cost_ul], 0.0)

## S2 -----------------------------------------------------------
function   S2()
  mds_symb = [:S2Mkt, :S2Def, :S2Life, :S2Health, :S2NonLife]
  mds = Array(S2Module, 0)
  mds_scr = zeros(Float64, length(mds_symb), 2)
  balance = DataFrame()
  corr = zeros(Float64, length(mds_symb), length(mds_symb))
  bscr = zeros(Float64, 2)
  adj_tp = 0.0
  adj_dt = 0.0
  op = S2Op(zeros(Float64, 5)...)
  scr = 0.0
  return(S2(mds_symb, mds, mds_scr, balance, corr,
            bscr, adj_tp, adj_dt, op, scr))
end

function  S2(p::ProjParam,
             ds2_mkt_int::Dict{Symbol, Any},
             ds2_mkt_eq::Dict{Symbol, Any},
             eq2type::Dict{Symbol, Symbol},
             ds2_mkt::Dict{Symbol, Any},
             ds2_op::Dict{Symbol, Float64},
             ds2::Dict{Symbol, Any})
  s2 = S2()
  s2.corr = ds2[:corr]
  s2.balance = s2bal(p)
  merge!(ds2_op,
         [:tp => s2.balance[1, :tpg] + s2.balance[1, :bonus]])
  s2.op = S2Op(ds2_op)

  push!(s2.mds, S2Mkt(p, s2.balance,
                      ds2_mkt_int,
                      ds2_mkt_eq, eq2type,
                      ds2_mkt))
  push!(s2.mds, S2Def(p, s2.balance))
  push!(s2.mds, S2Life(p, s2.balance))
  push!(s2.mds, S2Health(p, s2.balance))
  push!(s2.mds, S2NonLife(p, s2.balance))
  bscr!(s2)
  scr!(s2.op, s2.bscr[GROSS])
  scr!(s2)
  return s2
end


## other functions ==============================================
## S2 balance sheet (unshocked)
function s2bal(p::ProjParam)
  invs = InvPort(p.t_0, p.T, p.cap_mkt, p.invs_par...)
  proj = Projection(p.proj_par..., p.cap_mkt, invs,
                    p.l_ins, p.l_other, p.dyn)
  return hcat(proj.val_0, DataFrame(scen = :be))
end

## S2 balance sheet (shocked)
function s2bal(p::ProjParam,
               md::S2Module,
               shock!::Any,
               scen::Symbol)
  cpm = deepcopy(p.cap_mkt)
  l_ins = deepcopy(p.l_ins)
  if shock! == nothing
    invs = InvPort(p.t_0, p.T, cpm, p.invs_par...)
  else
    if md.shock_object == :CapMkt shock!(cpm, md) end
    invs = InvPort(p.t_0, p.T, cpm, p.invs_par...)
    if md.shock_object == :InvPort shock!(invs, md) end
    if md.shock_object == :LiabIns shock!(l_ins, md) end
    if md.shock_object == :InvPort_LiabIns
      shock!(cpm, l_ins, md)
    end
  end
  proj = Projection(p.proj_par..., cpm, invs,
                    l_ins, p.l_other, p.dyn)
  return hcat(proj.val_0, DataFrame(scen = scen))
end

## basic own funds
bof(md::S2Module, scen::Symbol) =
  md.balance[md.balance[:scen] .== scen, :invest][1,1] -
  md.balance[md.balance[:scen] .== scen, :tpg][1,1] -
  md.balance[md.balance[:scen] .== scen, :l_other][1,1] -
  md.balance[md.balance[:scen] .== scen, :bonus][1,1]

## future discretionary benefits
fdb(md::S2Module, scen::Symbol) =
  md.balance[md.balance[:scen] .== scen, :bonus][1,1]

## S2MktInt -----------------------------------------------------
function scr!(s2_mkt_int::S2MktInt)
  scr_net =
    bof(s2_mkt_int, :be) .-
  float64([bof(s2_mkt_int, sm) for sm in s2_mkt_int.shock_type])
  scr_gross =
    scr_net .- fdb(s2_mkt_int, :be) +
    float64([fdb(s2_mkt_int, sm) for sm in s2_mkt_int.shock_type])

  i_up = findin(s2_mkt_int.shock_type, [:spot_up])[1]
  i_down = findin(s2_mkt_int.shock_type, [:spot_down])[1]

  s2_mkt_int.scen_up = scr_net[i_up] >= scr_net[i_down]
  s2_mkt_int.scr[NET] = maximum([0.0, scr_net])
  s2_mkt_int.scr[GROSS] =
    max(0.0,
        s2_mkt_int.scen_up ? scr_gross[i_up] : scr_gross[i_down])
end

function mktintshock!(cap_mkt::CapMkt,
                      s2_mkt_int,
                      int_type::Symbol)

  len = min(length(cap_mkt.rfr.x),
            length(s2_mkt_int.shock[:spot_up]),
            length(s2_mkt_int.shock[:spot_down]))
  spot = forw2spot(cap_mkt.rfr.x[1:len])


  if int_type == :spot_down
    forw =
      spot2forw(spot .*
                (1 .+ s2_mkt_int.shock[:spot_down][1:len]))
  else
    forw =
      spot2forw(spot .+
                max(spot .* s2_mkt_int.shock[:spot_up][1:len],
                    s2_mkt_int.spot_up_abs_min))
  end
  cap_mkt.rfr.x = deepcopy(forw)
end

## S2MktEq ------------------------------------------------------
function scr!(s2_mkt_eq::S2MktEq)
  scr_net =
    bof(s2_mkt_eq, :be) .-
  float64([bof(s2_mkt_eq, sm) for sm in s2_mkt_eq.shock_type])
  scr_gross =
    scr_net .- fdb(s2_mkt_eq, :be) +
    float64([fdb(s2_mkt_eq, sm) for sm in s2_mkt_eq.shock_type])

  s2_mkt_eq.scr[NET] =
    sqrt(scr_net' * s2_mkt_eq.corr * scr_net)[1]
  s2_mkt_eq.scr[GROSS] =
    sqrt(scr_gross' * s2_mkt_eq.corr * scr_gross)[1]
end

function mkteqshock!(invs::InvPort, s2mkt_eq, eq_type::Symbol)
  invs.mv_0 -= invs.igs[:IGStock].mv_0
  invs.igs[:IGStock].mv_0 = 0.0
  for invest in invs.igs[:IGStock].investments
    if s2mkt_eq.eq2type[invest.name] == eq_type
      invest.mv_0 *= (1 - s2mkt_eq.shock[eq_type])
    end
    invs.igs[:IGStock].mv_0 += invest.mv_0
  end
  invs.mv_0 += invs.igs[:IGStock].mv_0
end

## S2Mkt --------------------------------------------------------
function scr!(s2_mkt::S2Mkt)
  scr_net = zeros(Float64, length(s2_mkt.mds))
  scr_gross = zeros(Float64, length(s2_mkt.mds))

  scen_up = false
  for i = 1:length(s2_mkt.mds)
    scr_net[i] = s2_mkt.mds[i].scr[NET]
    scr_gross[i] = s2_mkt.mds[i].scr[GROSS]
    if :scen_up in names(s2_mkt.mds[i])
      scen_up = s2_mkt.mds[i].scen_up
    end
  end
  corr = (scen_up ? s2_mkt.corr_up : s2_mkt.corr_down)
  s2_mkt.scr[NET] = sqrt(scr_net' * corr * scr_net)[1]
  s2_mkt.scr[GROSS] = sqrt(scr_gross' * corr * scr_gross)[1]
end

## S2Op ---------------------------------------------------------
function scr!(s2op::S2Op, bscr)
  scr_op_prem =
    0.04 *
    (s2op.prem_earned +
       max(0, 1.2 * (s2op.prem_earned - s2op.prem_earned_prev)))
  scr_op_tp = 0.0045 * max(0, s2op.tp)
  s2op.scr =
    min(0.3 * bscr,
        max(scr_op_prem, scr_op_tp)) + 0.25 * s2op.cost_ul
end

## S2 -----------------------------------------------------------
function bscr!(s2::S2)
  scr_net = zeros(Float64, length(s2.mds))
  scr_gross = zeros(Float64, length(s2.mds))
  for i = 1:length(s2.mds)
    scr_net[i] = s2.mds[i].scr[NET]
    scr_gross[i] = s2.mds[i].scr[GROSS]
  end
  s2.bscr[NET] =  sqrt(scr_net' * s2.corr * scr_net)[1]
  s2.bscr[GROSS] = sqrt(scr_gross' * s2.corr * scr_gross)[1]
end

function scr!(s2::S2)
  s2.adj_dt = 0.0 ## fixme: deferred tax not implemented
  s2.adj_tp =
    max(0.0, min(s2.bscr[GROSS] - s2.bscr[NET], fdb(s2, :be)))
  s2.scr = s2.bscr[GROSS] - s2.adj_tp - s2.adj_dt + s2.op.scr
end
