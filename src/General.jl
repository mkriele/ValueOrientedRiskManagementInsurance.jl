module ValueOrientedRiskManagementInsurance

using Distributions
using DataFrames, DataArrays

# import Base.show, Base.isequal
# import Base.merge!
import Distributions.rand

## Types

# Simple economic capital model
include("ECModel/ECModel__Types.jl")
include("ECModel/ECModel_Functions.jl")




end # module
