abstract type ZBasis end

function _indices end

"""
    basissize(basis)

Return number of basis vectors in the `basis`.

See:
- [`zbasis`](@ref)
- [`symmetrized_basis`](@ref)
"""
function basissize end

"""
    FullZBasis(N)

Represents the states of a system of N.

See also: [`zbasis`](@ref)
"""
struct FullZBasis <: ZBasis
    N::Int
    FullZBasis(N) = N > 0 ? new(N) : throw(ArgumentError("N=$N needs to be a positive integer!"))
end

_indices(fzb::FullZBasis) = 1:2^fzb.N

basissize(fzb::FullZBasis) = 2^fzb.N

"""
    ZBlockBasis(N, k)

Represents the states of a system of N spins whith k |↑⟩ (magnetization = (k-N)/2).

See also: [`zbasis`](@ref)
"""
struct ZBlockBasis <: ZBasis
    N::Int
    k::Int
    function ZBlockBasis(N, k)
        N > 0 || throw(ArgumentError("N=$N needs to be a positive integer!"))
        (0 <= k) && (k <= N) || throw(ArgumentError("k=$k needs to be between 0 and $N."))
        return new(N,k)
    end
end

function _indices(zbb::ZBlockBasis)
    N, k = zbb.N, zbb.k
    if k == 0 || k == N
        [2^k]
    end
    inds = ones(Int, binomial(N, k))# initialize with 1 because "index = binary representation + 1"
    _zblock_inds!(inds, N, k)
    return inds
end

basissize(zbb::ZBlockBasis) = binomial(zbb.N, zbb.k)

"""
    zbasis(N[, k])

Represent a full z-basis for N spins. If k is provided, this represents only the block
with k |↑⟩ (so magnetization of (k-N)/2).
"""
zbasis(N) = FullZBasis(N)
zbasis(N, k) = ZBlockBasis(N, k)

## build indices in-place recursively
## Note: Indices computed are from 0 to 2^N -1 -> need to +1
## Implementation based on recursive factorization:
## - {|N, k⟩} = |↑⟩ ⊗ {|N-1,k⟩} ∪ |↓⟩ ⊗ {|N-1,k-1⟩}
## - {|N, N⟩} = {|2^N -1⟩}
## - {|N, 0⟩} = {|0⟩}
## - |1⟩ ⊗ |s⟩ = |1*2^N + s⟩ (|s⟩ has N spins)
function _zblock_inds!(states, N, k)
    if k == 0
        states
    elseif N == k
        states .+= 2^N - 1
    else
        front = binomial(N-1, k)
        @views states[front+1:end] .+= 2^(N-1)
        _zblock_inds!(view(states, 1:front), N-1, k)
        _zblock_inds!(view(states, front+1:length(states)), N-1, k-1)
        states
    end
end

"""
    SymmetrizedBasis

Not intended for direct use. See [`symmetrized_basis`](@ref).
"""
struct SymmetrizedBasis
    basis::ZBasis
    symmetries::Vector{AbstractSymmetry}
    sectors::Vector{Int}
end

"""
    symmetrized_basis(N[, k], symmetry, sector, more...)
    symmetrized_basis(zbasis, symmetry, sector, more...)

Construct a basis in the specified symmetry sectors. Any number of symmetries may be specified.

Either provide number of spins (and optionally `k` block) or a [`zbasis`](@ref).
"""
function symmetrized_basis(N::Int, symmetry::AbstractSymmetry, sector::Int, more...)
    symmetrized_basis(zbasis(N), symmetry, sector, more...)
end

function symmetrized_basis(N::Int, k::Int, symmetry::AbstractSymmetry, sector::Int, more...)
    symmetrized_basis(zbasis(N, k), symmetry, sector, more...)
end

function symmetrized_basis(zbasis::ZBasis, symmetry::AbstractSymmetry, sector::Int, more...)
    mod(length(more), 2) == 0 || ArgumentError("Odd number of arguments. Please provide 1 sector for each symmetry.")
    SymmetrizedBasis(zbasis, [symmetry, more[1:2:end]...], [sector, more[2:2:end]...])
end

basissize(basis::SymmetrizedBasis) = length(_phase_factors(_indices(basis.basis), basis.symmetries, basis.sectors))

symmetrize_state(state, args...) = symmetrize_state(state, symmetrized_basis(args...))

"""
    symmetrize_state(state, basis)
    symmetrize_state(state, args...)

Symmetrize the given `state` into the symmetric sector specified by the [`symmetrized_basis`](@ref).

Alternatively, provide everything needed to construct the [`symmetrized_basis`](@ref) and
will be constructed internally.
"""
function symmetrize_state(state, basis::SymmetrizedBasis)
    if length(state) != 2^basis.basis.N
        throw(ArgumentError("""State has wrong size.
            Expected $(2^basis.basis.N), got: $(length(state))"""))
    end

    inds = _indices(basis.basis)
    factors = _phase_factors(inds, basis.symmetries, basis.sectors)
    result = Vector{ComplexF64}(undef, length(factors))

    for (i, component) in enumerate(factors)
        tmp = ComplexF64(0)
        for (index, phase) in component
            tmp += exp(im*2π*phase)*state[index]
        end
        result[i] = tmp / √(length(component))
    end
    return result
end

symmetrize_operator(operator, args...) = symmetrize_operator(operator, symmetrized_basis(args...))

"""
    symmetrize_operator(operator, basis)
    symmetrize_operator(operator, args...)

Symmetrize the given `operator` into the symmetric sector specified by the [`symmetrized_basis`](@ref).

Alternatively, provide everything needed to construct the [`symmetrized_basis`](@ref) and
will be constructed internally.
"""
function symmetrize_operator(operator, basis::SymmetrizedBasis)
    if length(operator) != 4^basis.basis.N || size(operator, 1) != size(operator, 2)
        throw(ArgumentError("""Operator has wrong size.
            Expected $(2^basis.basis.N)x$(2^basis.basis.N), got: $(size(operator))"""))
    end
    inds = _indices(basis.basis)
    factors = _phase_factors(inds, basis.symmetries, basis.sectors)
    result = Matrix{ComplexF64}(undef, length(factors), length(factors))

    for (j, component2) in enumerate(factors)
        for (i, component1) in enumerate(factors)
            tmp = ComplexF64(0)
            for (index2, phase2) in component2
                for (index1, phase1) in component1
                    tmp += exp(im*2π*(phase1 - phase2))*operator[index1, index2]
                end
            end
            result[i,j] = tmp / √(length(component1)) / √(length(component2))
        end
    end
    return result
end


function _compatible(d1, d2)
    common_keys = intersect(keys(d1), keys(d2))
    # can check for equality since phases are normalized to the interval [0,1)
    # in _phase_factors!
    return all(k -> d1[k] == d2[k], common_keys)
end

"""
    _phase_factors(inds, symms, sectors)

Compute the phase factors needed for a basis transformation.

# Returns
- `Vector[Dict{eltype(inds), Rational{Int}}]`: Each entry in the vector corresponds to a
basis vector in the new basis. The `Dict` holds the contributing basis vectors of the old
basis with its phase factor in units of 2π.
"""
function _phase_factors(inds, symms, sectors)
    # How this algorithm works:
    #
    # Setting:
    # Imagine all basis vectors as the vertices of a graph and the symmetries generate
    # (directed) edges between them. These edges carry a phase factor that's exp(i2π*k/L),
    # where k is the symmetry sector and L the cycle_size of the symmetry.
    # To find the basis transformation from the old basis to the new symmetrized basis,
    # we need to find the connected components. Each connected component corresponds to a
    # new basis vector. All vectors in a component contribute equally (by magnitude) and to
    # find the phase factors we just follow the edges.
    # This is complicated by the fact, that the graph is not acyclic! Thus we need to make
    # sure that the different symmetries are actually compatible in the component. That
    # means that every starting and ending at specified nodes needs to have the same phase.
    # Otherwise the whole component just vanishes under the combined symmetry.
    #
    # Implementation detail:
    # We use Rational numbers r instead of Complex64 φ to track the phases. Their relation is
    # φ = exp(i2π*r)
    #
    # Algorithm:
    # 1. Start at a basis vector
    # 2. Find connected vectors by applying all symmetries and checking for compatibility
    # of their phases
    # 3. Track the phases of all vectors seen so far (check for compatibility of phases)
    # 4. Continue generating connected vectors for all vectors until there are no new discoveries
    # 5. Repeat all steps for the next unseen basis vector
    DTYPE = Dict{eltype(inds), Rational{Int}}
    seen = Set{eltype(inds)}()
    output = DTYPE[]
    for i in inds
        i ∈ seen && continue # vector already taken care of

        out = Dict(i => 0//1)
        active = Dict(i => 0//1)
        lives = true
        oldcount = 0 # track no. of vectors found
        while oldcount != length(out)
            oldcount = length(out)

            new_active = DTYPE()
            for (k,s) in zip(sectors, symms)
                add = DTYPE(apply(s, i)=>mod(c+k//_order(s),1) for (i,c) in active)
                lives = lives && _compatible(add, new_active)
                merge!(new_active, add)
            end
            active = new_active

            lives = lives && _compatible(active, out)

            merge!(out, active)
        end
        union!(seen, keys(out))
        lives && push!(output, out)
    end
    output
end
