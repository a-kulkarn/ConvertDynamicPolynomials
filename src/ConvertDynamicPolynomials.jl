
module ConvertDynamicPolynomials

import AbstractAlgebra
import Nemo
import DynamicPolynomials

# aliases
DP = DynamicPolynomials

export coeftype
coeftype(::Type{DP.Polynomial{C, T}}) where {C, T} = T
coeftype(p::DP.Polynomial{C, T}) where {C, T} = T

Base.one(X::Vector{DP.PolyVar{true}}) = monomials(X,0)[1]


# Macro for the old algebraic solvers interface.
# Generates a DynamicPolynomial "ring".
export @Ring
macro Ring(args...)
    X = DP.PolyVar{true}[DP.PolyVar{true}(string(arg)) for arg in args]
    V = [buildpolyvar(DP.PolyVar{true}, args[i], X[i]) for i in 1:length(X)]
    push!(V, :(TMP = $X) )
    reduce((x,y) -> :($x; $y), V; init = :() )
end

function buildpolyvar(::Type{PV}, arg, var) where PV
    :($(esc(arg)) = $var)
end


GeneralDP = DP.Polynomial{B, S} where {B, S}

## Need new versions of the AA_* constructors for terms/monomials/etc...


export DPmonomials
DPmonomials = DP.monomials

export AAPolynomialRing, AAPolynomial

function extract_variables(P)
    mons = P.x # Array of variables of P
    return unique(vcat( [m.vars for m in mons]...))    
end

"""
    Given a coefficient ring, which is a subtype of NCRing, return the AbstractAlgebra
"""
function AAPolynomialRing(coeff_ring::T where T <: AbstractAlgebra.NCRing,
                          X::Array{DP.PolyVar{B},1} where B )
    
    return PolynomialRing(coeff_ring, [string(x) for x in X])
end

function AAPolynomialRing(P::T where T <:GeneralDP )
    if DP.maxdegree(P)==0
        error("Cannot construct AbstractAlgebra parent from degree 0 polynomial")
    elseif !(typeof(P.a[1]) <: AbstractAlgebra.NCRing) && typeof(P.a[1]) != Int64
        error("Coefficient type has no canonical parent")
    end

    if typeof(P.a[1]) == Int64
        R = Nemo.FlintIntegerRing()
    else
        R = parent(a[1])
    end
    
    vars = extract_variables(P)

    return AbstractAlgebra.PolynomialRing(R, [string(v) for v in vars])
end

"""

    Converts a DynamicPolynomial to an AbstractAlgebra polynomial with parent `RX`. If no parent is
    provided, once is guessed using  "AAPolynomialRing".

    Default assignment of Dynamic Polynomial variables to the parent is greedy. However, one
    can provide a dictionary to make explicit the desired assignment.

"""

function AAPolynomial(P::T where T <:GeneralDP )
    future_parent, x = AAPolynomialRing(P)
    AAPolynomial(P, future_parent)
end


# TODO: type assert the var_assignment input at some point.
function AAPolynomial(P::T where T <:GeneralDP , RX; var_assignment=nothing)

    if var_assignment==nothing

        dvars = extract_variables(P)

        if size(dvars,1) > size(AbstractAlgebra.gens(RX),1)
            error("Input polynomial has more variables than adoptive parent ring")
        end

        varDic = Dict( zip(dvars,  AbstractAlgebra.gens(RX)[1:length(dvars)]) )
    else
        varDic = var_assignment
    end
    
    coeffs = [RX(c) for c in P.a]

    function coerce_monomial(m)
        e = m.z
        v = [ varDic[x] for x in m.vars]
        r = size(v,1)

        if r > 0
            return prod( v[i]^e[i] for i=1:r)
        else
            return RX(1)
        end
    end

    
    mons = [coerce_monomial(m) for m in P.x]
    r = size(mons,1)
    
    if iszero(P)
        return zero(RX)
    else
        return sum( mons[i]*coeffs[i] for i=1:r ) 
    end

end


## Allegedly a convenience function to evaluate polyonmials by P(X).
## Originally part of AlgebraicSolvers
## BUG: It allows you to evaluate at a vector that is too long...
# function (p::DP.Polynomial{B,T})(x::Vector) where {B,T}
#    r = zero(x[1]);
#    for m in p
#       t=m.α
#       for i in 1:length(m.x.z)
#       	 t*=x[i]^m.x.z[i]
#       end
#       r+=t
#    end
#    r
# end



    # if vars==nothing && parent==nothing
    #     if iszero(P) || !( typeof(P.a[1]) <: Hecke.NCRingElem )
    #         error("Cannot convert constant polynomial without an explicit parent for the coefficients.")
    #     end
    # end

    # if parent==nothing
    #     RX, RXvars = AAPolynomialRing(parent(P.a[1]), vars)
    # else
    #     vars = extract_variables(P):: T where T <: Array{DP.PolyVar{true},1}
    #     RX, RXvars = AAPolynomialRing(parent(P.a[1]), vars)
    # end

end # module
