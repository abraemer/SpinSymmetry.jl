# ToDo
- nothing concrete
- more symmetry operations? lattice symmetries maybe?

# Further ideas
- for interactive development: Get representants of the vectors in the symmetrized basis
- optimize speed/allocations? Hardly worth it. It's quite fast as it is.

# Concerning speed
Quite fast... While constructing the phases is fast, the actual symmetrization is quite horrible. The easy answer is to convert the phases to a SparseArray and simply multiply. This gives several orders of magnitude speedup and apparently does not influence the allocations much. It also may be cached to if needed to speed up symmetrization even more.

Before:
```julia
julia> @btime symmetrize_operator(m, 14, Flip(14),0)
  37.104 s (385096 allocations: 53.52 MiB)
```

After:
```julia
julia> m2 = @btime SpinSymmetry.symmetrize_operator2(m, symmetrized_basis(14, Flip(14),0))
  37.654 ms (385142 allocations: 55.91 MiB)

julia> m1 ≈ m2
true
```