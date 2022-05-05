# SpinSymmetry

[![codecov](https://codecov.io/gh/abraemer/SpinSymmetry.jl/branch/main/graph/badge.svg?token=XN6TT95A53)](https://codecov.io/gh/abraemer/SpinSymmetry.jl)
[![CI](https://github.com/abraemer/SpinSymmetry.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/abraemer/SpinSymmetry.jl/actions/workflows/ci.yml)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Transform your spin system to a symmetry sector of your choice!

SpinSymmetry.jl is a light-weight, reasonably fast package with no dependencies. Some commonly used symmetries are implemented, there's support for user-defined symmetries and every symmetry can be combined!

# Install
```julia
julia> Pkg.add("SpinSymmetry")
```

# Usage
Construct a `SymmetrizedBasis` and use it to transform your state vectors and operators:
```julia
julia> using SpinSymmetry, LinearAlgebra, BenchmarkTools

julia> basis = symmetrized_basis(zbasis(9, 4), Flip(9), 0, Shift(9), 0);

julia> state = normalize!(ones(2^9)); # initial state - all up in x direction

julia> @btime symm_state = symmetrize_state(state, basis)
637.317 μs (4245 allocations: 646.05 KiB)
14-element Vector{Float64}:
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994
 0.18749999999999994

julia> operator = kron([1 0; 0 -1], kron([1 0; 0 -1], I(2^7))); # operator Z ⊗ Z ⊗ 𝟙 ⊗ .. ⊗ 𝟙

julia> @btime symm_op = symmetrize_operator(operator, basis)
 1.096 ms (4248 allocations: 703.76 KiB)
14×14 Matrix{Float64}:
[...]
```

# Features
- Use `symmetrized_basis` to construct a collection of your symmetries. Provide as the first argument the system's size (and optionally magnetizaion block) and then follow with symmetry operations and their sector alternating.
- Apply the symmetries to a state or an operator using `symmetrize_state` and `symmetrize_operator`
- Find the size of the symmetry sector with `basissize`

The symmetry operations supported are:
- z-magnetization block (via `zbasis(N, k)`)
- Spin flip via `Flip(positions)` or `Flip(N)`
- Shift symmetry via `Shift(N, amount=1)`
- Swap/Exchange symmetry via `Swap(pos1, pos2)`
- Spatial reflection via `SpatialReflection(N)`

where `N` denotes the number of spins in the system and their positions should be given as a Julian index, i.e. in the range `1:N`.

To get the tranformation matrix to the symmetrized subspace just use `transformationmatrix(symmetrized_basis)`.

**Note:** The projection on a specific magnetization block is applied first. Thus if you have spin flip symmetry and restrict to a magnetization block, your symmetrized basis states look like "|↑..↑⟩ ± |↓..↓⟩". So in this case you effectively specified S_z^2 and parity.

## User-defined symmetries
It's also quite easy to define your own symmetry operations. 
Simply define a function `f` that maps one basis index to the next.
Note that these basis indices are a binary representation (range 0:2^N-1) of the spin basis where the first spin is represented by the *least* significant bit.
Then you can use `GenericSymmetry(f, L)` where `L` denotes the order of your symmtry.
The order is the smallest number `L` s. t. `f` applied `L` is the identity for all indices.

Suppose the spatial reflection would not be implemented. You could do it yourself by defining:
```julia
julia> reflection(N) = bits -> parse(Int, string(bits; base=2, pad=N)[end:-1:1]; base=2)
julia> SpatialReflection(N) = GenericSymmetry(reflection(N), 2)
```

# Implementation details
Imagine all basis vectors as the vertices of a graph and the symmetries generate
(directed) edges between them. These edges carry a phase factor that's `exp(i2π*k/L)`,
where `k` is the symmetry sector and `L` the order of the symmetry.
To find the basis transformation from the old basis to the new symmetrized basis,
we need to find the connected components. Each connected component corresponds to a
new basis vector. All vectors in a component contribute equally (by magnitude) and to
find the phase factors we just follow the edges.
This is complicated by the fact, that the graph is not acyclic! Thus we need to make
sure that the different symmetries are actually compatible in the component. That
means that every starting and ending at specified nodes needs to have the same phase.
Otherwise the whole component just vanishes under the combined symmetry.

Implementation detail:
We use Rational numbers `r` instead of Complex64 `φ` to track the phases. Their relation is
`φ = exp(i2π*r)`

Algorithm:
1. Start at a basis vector
2. Find connected vectors by applying all symmetries and checking for compatibility
of their phases
3. Track the phases of all vectors seen so far (check for compatibility of phases)
4. Continue generating connected vectors for all vectors until there are no new discoveries
5. Repeat all steps for the next unseen basis vector
