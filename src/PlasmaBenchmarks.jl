module PlasmaBenchmarks

using Plasma
using BSON: @save

include("two-stream.jl")
include("weibel.jl")
# include("1d1v.jl")
# include("3d3v.jl")

end