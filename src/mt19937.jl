# see also https://en.wikipedia.org/wiki/Mersenne_Twister
# see also cryptopals.com (2026) Set 3 for untemper and clonefromoutout rationale

module MT19937

import Base: rand

export MT19937State, initialize_state!, random_uint32!, rand_float64!, MersenneTwister
export untemper, clonefromoutput

const N = 624
const M = 397
const MATRIX_A = UInt32(0x9908b0df)
const UPPER_MASK = UInt32(0x80000000)
const LOWER_MASK = UInt32(0x7fffffff)

""" MT19937 internal state structure: an array of 624 UInt32 and an index. """
mutable struct MT19937State
    mt::Vector{UInt32}
    index::Int
    function MT19937State(a = zeros(UInt32, N), b = 0)
        new(a, b)
    end
end

""" Initialize the MT19937 state with a UUInt32 seed. """
function initialize_state!(state::MT19937State, seed::UInt32)
    state.mt[1] = seed
    for i in 2:N
        state.mt[i] = UInt32(1812433253) *
                      (state.mt[i-1] ⊻ (state.mt[i-1] >> 30)) +
                      UInt32(i-1)
    end
    state.index = 0
end

""" Twist the internal state of the MT19937 RNG. """
function mt_twist!(rng::MT19937State)
    for i in 1:N
        x = (rng.mt[i] & UPPER_MASK) |
            (rng.mt[mod1(i, N)] & LOWER_MASK)

        xA = x >> 1
        if (x & 1) != 0
            xA ⊻= MATRIX_A
        end

        rng.mt[i] = rng.mt[mod1(i + M, N)] ⊻ xA
    end
    rng.index = 0
end

""" Generate a random UInt32 using MT19937. """
function random_uint32!(rng::MT19937State)::UInt32
    if rng.index ≥ N
        mt_twist!(rng)
    end

    y = rng.mt[rng.index+1]
    rng.index += 1

    # Tempering (exact)
    y ⊻= y >> 11
    y ⊻= (y << 7) & UInt32(0x9D2C5680)
    y ⊻= (y << 15) & UInt32(0xEFC60000)
    y ⊻= y >> 18

    return y
end

""" Generate a random Float64 in [0,1) using MT19937. """
function rand_float64!(rng::MT19937State)::Float64
    a = random_uint32!(rng) >> 5
    b = random_uint32!(rng) >> 6
    return (Float64(a) * 67108864.0 + Float64(b)) / 9007199254740992.0
end

""" Mersenne Twister MT19937 RNG wrapper. """
struct MersenneTwister
    state::MT19937State
    function MersenneTwister(state::MT19937State)
        new(state)
    end
    function MersenneTwister(seed::UInt32)
        state = initialize_state!(MT19937State(zeros(UInt32, N), 0), seed)
        new(state)
    end
end
Base.rand(mt::MersenneTwister) = rand_float64!(mt.state)

""" Invert the tempering function to recover one portion of internal state from an output. """
function untemper(y::UInt32)::UInt32
    # Invert y ^= y >> shift
    function unshift_right_xor(y::UInt32, shift::Int)
        x = y
        for _ in 1:ceil(Int, 32/shift)
            x = y ⊻ (x >> shift)
        end
        return x
    end
    # Invert y ^= (y << shift) & mask
    function unshift_left_xor_mask(y::UInt32, shift::Int, mask::UInt32)
        x = y
        for _ in 1:ceil(Int, 32/shift)
            x = y ⊻ ((x << shift) & mask)
        end
        return x
    end

    x = y
    x = unshift_right_xor(x, 18)
    x = unshift_left_xor_mask(x, 15, 0xEFC60000)
    x = unshift_left_xor_mask(x, 7, 0x9D2C5680)
    x = unshift_right_xor(x, 11)
    return x
end

""" Clone the internal state from 624 outputs, as per Cryptopals exercise 22. """
function clonefromoutput(outputs::Vector{UInt32})
    @assert length(outputs) == N
    mt = Vector{UInt32}(undef, N)
    for i in 1:N
        mt[i] = untemper(outputs[i])
    end
    return MT19937State(mt, N)
end


end # module MT19937
