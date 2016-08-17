using ValueOrientedRiskManagementInsurance
using DataFrames
using Base.Test

println("start ECModel2.jl...")

lob_strings = ["Fire", "Liability", "Theft"]
df_claims = Vector(length(lob_strings))
claims = Vector(length(lob_strings))
mack = Vector(length(lob_strings))
β = Vector(length(lob_strings))
for i = 1:length(df_claims)
  df_claims[i] =
    readtable("test/ECModel_Input_" * lob_strings[i] * ".csv",
              header = false)
  mack[i] = Mack(df_claims[i])
  β[i] = mack[i].futureclaims / sum(mack[i].futureclaims)
end

function logpar2statpar(m_s::Vector)
  Real[exp(m_s[1] + m_s[2]^2 / 2),
       √(exp(m_s[2]^2-1)) * exp(m_s[1] + m_s[2]^2 / 2)]
end

function statpar2logpar(μ_σ::Vector)
  vc = μ_σ[2]/μ_σ[1]
  Real[log(μ_σ[1]) - 0.5 * (1+log(vc^2)), √(1+log(vc^2))]
end

logpar2statpar([0,1])
x = logpar2statpar([0,1])


statpar2logpar(x)

statpar2logpar(x[1],x[2])
#
# mean, sde -> logmean, logsigma
# R[i] = LogNormal(log_mean[i], log_sigma[i])
# v⋅mack[i].β * R[i]
#
# v: stochastic discont
# β
# R[i] = stochastic total undiscounted future paiments,
#

println("...end ECModel2.jl")
