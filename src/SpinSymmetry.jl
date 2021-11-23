module SpinSymmetry

using SparseArrays: sparse

export Flip, Shift, Swap, SpatialReflection, GenericSymmetry
export zbasis, FullZBasis, ZBlockBasis, basissize
export SymmetrizedBasis, symmetrized_basis, symmetrize_state, symmetrize_operator, transformationmatrix

include("abstract.jl")
include("basis.jl")
include("symmetries.jl")

end
