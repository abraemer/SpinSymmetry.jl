module SpinSymmetry

import SparseArrays

export Flip, Shift, Swap, SpatialReflection, GenericSymmetry
export zbasis, FullZBasis, ZBlockBasis, basissize
export SymmetrizedBasis, symmetrized_basis, symmetrize_state, symmetrize_operator

include("abstract.jl")
include("basis.jl")
include("symmetries.jl")

end
