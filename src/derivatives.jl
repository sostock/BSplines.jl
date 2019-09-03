"""
    AllDerivatives{N}

A singleton type that represents all `m`-th derivatives where `0 ≤ m < N`.
"""
struct AllDerivatives{N}
    function AllDerivatives{N}() where N
        N isa Int || throw(ArgumentError("Type parameter N in AllDerivatives{N} must be an Int"))
        N ≥ 1 || throw(DomainError(N, "Type parameter N in AllDerivatives{N} must be positive"))
        new()
    end
end

"""
    AllDerivatives(N)

A shortcut for `AllDerivatives{N}()`, representing all `m`-th derivatives where `0 ≤ m < N`.
"""
AllDerivatives(N::Integer) = AllDerivatives{Int(N)}()

"""
    Derivative{N}

A singleton type that represents the `N`-th derivative.
"""
struct Derivative{N}
    function Derivative{N}() where N
        N isa Int || throw(ArgumentError("Type parameter N in Derivative{N} must be an Int"))
        N ≥ 0 || throw(DomainError(N, "Type parameter N in Derivative{N} must be positive"))
        new()
    end
end

"""
    Derivative(N)

A shortcut for `Derivative{N}()`, representing the `N`-th derivative.
"""
Derivative(N::Integer) = Derivative{Int(N)}()

"""
    NoDerivative

An alias for `Derivative{0}`. (not exported)
"""
const NoDerivative = Derivative{0}

"""
    NoDerivUnion

An alias for `Union{Derivative{0}, AllDerivatives{1}}`. (not exported)
"""
const NoDerivUnion = Union{Derivative{0}, AllDerivatives{1}}

# AllDerivatives/Derivative act as scalars for broadcasting
Base.Broadcast.broadcastable(x::AllDerivatives) = Ref(x)
Base.Broadcast.broadcastable(x::Derivative) = Ref(x)
