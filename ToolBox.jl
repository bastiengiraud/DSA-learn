module ToolBox


# specify current path
cd("C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code")

# activate path and show active packages
using Pkg
Pkg.activate(".")
Pkg.status()

# include support scripts
include("functions/method.jl")
include("functions/contingency.jl")
include("functions/ssa_module.jl")
include("functions/dynamics.jl")
include("functions/directed_walk.jl")
include("functions/acpfcorrect.jl")
include("functions/obbt_lu.jl")
include("functions/polytope.jl")
include("functions/support.jl")
include("functions/write_dfs.jl")


end
