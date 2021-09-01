# SpinSymmetry

[![codecov](https://codecov.io/gh/abraemer/SpinSymmetry.jl/branch/main/graph/badge.svg?token=XN6TT95A53)](https://codecov.io/gh/abraemer/SpinSymmetry.jl)
[![CI](https://github.com/abraemer/SpinSymmetry.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/abraemer/SpinSymmetry.jl/actions/workflows/ci.yml)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Transform your spin system to a symmetry sector of your choice!

# Install
```julia
julia> Pkg.add(url="https://github.com/abraemer/SpinSymmetry.jl")
```

# Usage
Construct a `SymmetrizedBasis` and use it to transform your state vectors and operators:
```julia
julia> using SpinSymmetry, LinearAlgebra, BenchmarkTools

julia> basis = symmetrized_basis(zbasis(9, 4), Flip(9), 0, Shift(9), 0);

julia> state = normalize!(ones(2^9)); # initial state - all up in x direction

julia> @btime symm_state = symmetrize_state(state, basis)
 765.253 μs (7781 allocations: 727.29 KiB)
14-element Vector{ComplexF64}:
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im
 0.18750000000000006 + 0.0im

julia> operator = kron([1 0; 0 -1], kron([1 0; 0 -1], I(2^7))); # operator Z ⊗ Z ⊗ 𝟙 ⊗ .. ⊗ 𝟙

julia> @btime symm_op = symmetrize_operator(operator, basis)
 41.828 ms (922891 allocations: 29.66 MiB)
14×14 Matrix{ComplexF64}:
[...]
```

# Features
The symmetry operations supported are:
- z-magnetization block (via `zbasis(N, k)`)
- Spin flip via `Flip(positions)` or `Flip(N)`
- Shift symmetry via `Shift(N, amount=1)`
- Swap/Exchange symmetry via `Swap(pos1, pos2)`

where `N` denotes the number of spins in the system and their positions should be given as a Julian index, i.e. in the range `1:N`.

**Note:** The projection on a specific magnetization block is applied first. Thus if you have spin flip symmetry and restrict to a magnetization block, your symmetrized basis states look like "|↑..↑⟩ ± |↓..↑⟩".

It's also quite easy to define your own symmetry operations. 
Simply define a function `f` that maps one basis index to the next.
Note that these basis indices are a binary representation of the spin basis where the first spin is represented by the *least* significant bit.
Then you can use `GenericSymmetry(f, L)` where `L` denotes the order of your symmtry.
The order is the smallest number `L` s. t. `f` applied `L` is the identity for all indices.

## Other functions
- Use `symmetrized_basis` to construct a collection of your symmetries
- If you only need them once, you can pass them to `symmetrize_state` and `symmetrize_operator` directly (instead of the `SymmetrizedBasis` object).

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