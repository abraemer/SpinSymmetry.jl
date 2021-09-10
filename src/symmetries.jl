"""
    Shift(N, [k=1])

Shift operator. Optionally specify the distance to shift (default: 1).

# Fields
- `N::Int`: total number of spins
- `amount::Int`: shift by how many places
"""
struct Shift <: AbstractSymmetry
    N::Int
    amount::Int
    Shift(N, k=1) = new(N, k)
end

_order(s::Shift) = div(s.N, gcd(s.N, s.amount))

_mask(k) = (2^k)-1 # mask out the k lowest bits
_roll_bits(N, i, amount) = ((i & _mask(amount)) << (N-amount)) | (i >> amount)
(s::Shift)(index) = _roll_bits(s.N, index, s.amount)


"""
    Flip(positions)
    Flip(system_size)

Spin Flip operator on N spins. Either provide the specific locations to flip or the system's
size to flip them all.

# Fields
- `mask::Int`: encodes the positions to flip in binary
"""
struct Flip <: AbstractSymmetry
    mask::Int
    """
        Flip(positions)
        Flip(system_size)

    Flip the spins at the positions or all if the system's size is given.
    """
    Flip(positions) = new(_compute_mask(positions))
end
Flip(system_size::Int) = Flip(1:system_size)
_compute_mask(inds) = sum(1 .<< (inds .- 1))

_order(::Flip) = 2

(f::Flip)(index)  = xor(index, f.mask)


"""
    Swap(pos1, pos2)

Swaps the 2 specified spins.

# Fields
- `ind1::Int`: Position of the first spin to flip
- `ind2::Int`: Position of the second spin to flip
Positions should be given in the range 1:N.
"""
struct Swap <: AbstractSymmetry
    ind1::Int
    ind2::Int
    function Swap(pos1, pos2)
        pos1 == pos2 && @warn "Swap with identical positions is nonsensical and won't work properly!"
        new(pos1, pos2)
    end
end

_order(::Swap) = 2

function _swap_bits(x, ind1, ind2)
    # compare bits, swap if unequal
    sw = ((x >> ind1) & 1) ⊻ ((x >> ind2) & 1) # xor == 1 if unequal
    return xor(x, (sw << ind1) | (sw << ind2)) # xor(bit, 1) -> flips bit
end

# convert spin positions to bitstring position
(s::Swap)(index) = _swap_bits(index, s.ind1-1, s.ind2-1)

"""
    SpatialReflection(N)

Reflects the whole chain in space.

# Fields
- `N::Int`: Number of spins
"""
struct SpatialReflection <: AbstractSymmetry
    N::Int
    function SpatialReflection(N)
        N == 1 && @warn "SpatialReflection with N=1 is nonsensical and won't work properly!"
        new(N)
    end
end

_order(::SpatialReflection) = 2

# approx 10times faster than the simpler
# parse(Int, string(x; base=2, pad=N)[N:-1:1]; base=2)
# Also note that the former does not account for cases where there are actually more than N
# spins
function _reflect_bits(N, bits)
    for i in 0:div(N,2)-1
        bits = _swap_bits(bits, i, N-1-i)
    end
    return bits
end

(sr::SpatialReflection)(index) = _reflect_bits(sr.N, index)

"""
    GenericSymmetry(f, L)

Allows a user-defined symmetry. Just provide a function to transform a binary basis index
and the order of the symmetry.

# Fields
- `f::Function`: transformation of the indices given in binary representation (0:2^N-1)
- `L::Int`: order of the symmetry, i.e. the smallest positive number such that f^L = identity
"""
struct GenericSymmetry <: AbstractSymmetry
    f::Function
    L::Int
end

_order(gs::GenericSymmetry) = gs.L

(gs::GenericSymmetry)(index) = gs.f(index)
