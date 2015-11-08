module ValueOrientedRiskManagementInsurance

export es

using Distributions
using DataFrames, DataArrays

# import Base.show, Base.isequal
# import Base.merge!
import Distributions.rand

## General functions --------------------------------------------
## Expected shortfall
es(x::Vector,α) =
  mean(sort(x, rev = true)[1:ceil(Integer, (1 - α) * length(x))])

# Simplfied Swiss Solvency Test----------------------------------
include("SST/SST__Types.jl")
include("SST/SST_Functions.jl")

# Simplified Life insurer ---------------------------------------
include("Life/Life__Types.jl")
include("Life/Life_Constructors.jl")
include("Life/Life_Functions.jl")

# Simplified Solvency 2 Life ------------------------------------
include("S2Life/S2Life__Types.jl")
include("S2Life/S2Life_Constructors.jl")
include("S2Life/S2Life_Functions.jl")

# Simplified Solvency 2 Non-Life --------------------------------
include("S2NonLife/S2NonLife__Types.jl")
include("S2NonLife/S2NonLife_Functions.jl")

# Simple economic capital model ---------------------------------
include("ECModel/ECModel__Types.jl")
include("ECModel/ECModel_Functions.jl")

end # module
