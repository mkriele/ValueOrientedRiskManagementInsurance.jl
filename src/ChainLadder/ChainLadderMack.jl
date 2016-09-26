using DataFrames

export Mack, claims2cum, cum2claims, cum2futureclaims, replacena!
"""
A reserve triangle with chain ladder information
"""
type Mack
  "Length of history"
  I::Int
  "Cummulative claims (upper triangle)"
  cum::Matrix{Real}
  "Time vector of development factors"
  f::Vector{Real}
  "Time vector of future claims payments"
  futureclaims::Vector{Real}
  "Vector of reseves"
  res::Vector{Real}
  """
  Total reserves: ``\mathit{tot\_res} =\sum_{i=1}^I \mathit{res}_i``
  """
  tot_res::Real
  "Vector of mean square errors"
  mse::Vector{Real}
  "Total mean square error"
  tot_mse::Real
end


# function Mack(df_triang::DataFrame; cum::Bool = false)
# Mack(convert(DataArray{Real, 2}, df_triang); cum = cum)
# end
#
# function Mack(da_triang::DataArray{Real, 2}; cum::Bool = false)
#   Mack(convert(Array{Real,2}, da_triang, 0); cum = cum)
# end

function Mack(triang::Array{Real,2}; cum::Bool = false)
  if cum
    C = deepcopy(triang)
  else
    C = claims2cum(triang)
  end
  # Number of accident / development years
  I = size(C)[1]
  # Development factors
  f = [sum(C[1:(I-k),k+1])/sum(C[1:(I-k),k]) for k in 1:(I-1)]
  # Ultimate claim amounts
  for k  in 2:I
    for l in (I+1-k+1):I
      C[k,l] = C[k,I+1-k] * prod(f[(I+1-k):(l-1)])
    end
  end
  futclaims = cum2futureclaims(C)
  # Reserves
  R = [C[k,I]-C[k,I+1-k] for k in 1:I]
  σ² =
    [1/(I-k-1) *
      C[1:(I-k),k] ⋅ (C[1:(I-k),k+1] ./ C[1:(I-k),k]-f[k]).^2
    for k in 1:(I-1)]
  σ²[I-1] = min(σ²[I-2]^2 / σ²[I-3], min(σ²[I-3],σ²[I-2]))
  # mean square error
  mse = zeros(Real, I)
  for i in 2:I
    mse[i] = 0
    for k in (I+1-i):(I-1)
      mse[i] += σ²[k]/f[k]^2 * (1/C[i,k] + 1/ sum(C[1:I-k,k]))
    end
    mse[i] *= C[i,I]^2
  end
  # total mean square error
  mse_total = 0.0
  for i in 2:I
    tmp = 0.0
    for k in collect(I+1-i : I-1)
      tmp += (2σ²[k] / f[k]^2) / sum(C[1 : I-k, k])
    end
    tmp *= C[i,I] * sum(C[i+1 : I, I])
    tmp += mse[i]
    mse_total += tmp
  end
  return Mack(I, C, f, futclaims, R, sum(R), mse, mse_total)
end

"""
replacena!(df::DataFrame, replacement::Any)

Replaces all 'NA'-values with 'replacement'
"""
function replacena!(df::DataFrame, replacement::Any)
  nrows, ncols = size(df)
  for j = 1:ncols; for i = 1:nrows
    if isna(df[i,j]); df[i,j] = replacement; end
    end
  end
end

"""
upperleft(quadrat_mat::Array{Real,2})

Sets all values in the strict lower right triangle to 0.
"""
function upperleft{T<:Real}(quadrat_mat::Array{T,2})
  c = deepcopy(quadrat_mat)
  n = size(c,1)
  if n ≠ size(c,2)
    error("known: Matrix is not quadratic")
  end
  for j = 1:n; for i = 1:n
    if i > n-j+1; c[i,j] = 0 end
  end; end
  c
end

"""
lowerright(quadrat_mat::Array{Real,2})

Sets all values in the upperleft triangle (incl. diagonal) to 0.
"""
function lowerright{T<:Real}(quadrat_mat::Array{T,2})
  quadrat_mat - upperleft(quadrat_mat)
end

"""
    claims2cum(c::Array{Real,2})

    Convert a quadratic upper left triangular reserve matrix to
    the correponding cumulative triangular reserve matrix.
    The (strictly) lower triangular part is set to 0
"""
function claims2cum(c::Array{Real,2})
  if size(c,1) ≠ size(c,2)
    error("claims2cum: Matrix is not quadratic")
  end
  cum = zeros(Real,size(c))
  #for i in 1:size(c, 1)
    #cum[i,1:size(c, 1)-i+1] = cumsum(c[i,1:size(c, 1)-i+1], 2)
  #end
  upperleft(cumsum(c, 2))
end

"""
    cum2claims(c::Array{Real,2})

    Convert a quadratic cumulative claims matrix into a
    claims matrix.
"""
function cum2claims(c::Array{Real,2})
  I = size(c,2)
  if size(c,1) ≠ I
    error("cum2claims: Matrix is not quadratic")
  end
  claims = zeros(Real,size(c))
  for i in 1:I
    claims[i,1] = c[i,1]
    for k in 2:I
     claims[i,k] = c[i,k]-c[i,k-1]
    end
  end
  claims
end

"""
    cum2futclaims(c::Array{Real,2})

    Calculate the future claims vector from the quadratic,
    (strictly) lower triangle of a projected cumulative claims
    matrix
"""
function cum2futureclaims(c::Array{Real,2})
  I = size(c,2)
  if size(c,1) ≠ I
    error("cumclaims: Matrix is not quadratic")
  end
  claims = cum2claims(c)
  futureclaims = zeros(Real,I-1)
  for k in 2:I
    futureclaims[k-1] = sum(claims[I-k+2:I,k])
  end
  futureclaims
end
