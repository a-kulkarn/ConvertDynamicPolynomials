
using AbstractAlgebra
using ConvertDynamicPolynomials

X = @Ring x1 x2 x3
n = length(X)-1


P = x1^2+1
display(P)
println(typeof(AAPolynomial(P)))
display(AAPolynomial(P))
println()

P = x1^0 + x2^0
println(typeof(AAPolynomial(P)))
