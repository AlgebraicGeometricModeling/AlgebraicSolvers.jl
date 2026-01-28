export  as_polynomial, convert_coeff, as_monomial, convert_DP

function as_monomial(p)
    DynamicPolynomials.monomials(p)[end]
end


function as_monomial(p, X)
    e = first(collect(exponent_vectors(p)))
    return prod(X[k]^e[k] for k in 1:length(e))
end

function as_polynomial(p, X, C::Type=Rational{BigInt})

    Exp   = collect(AbstractAlgebra.exponent_vectors(p))
    Coeff = [C(c) for c in AbstractAlgebra.coefficients(p)]
    
    Mon = [prod(X[k]^e[k] for k in 1:length(e)) for e in Exp]

    pol = sum(Coeff[i]*Mon[i] for i in 1:length(Coeff))

    return pol
end

function as_polynomial(p::DynamicPolynomials.Polynomial, R)

    Exp   = DynamicPolynomials.exponents.(DynamicPolynomials.monomials(p))
    Coeff = [c for c in DynamicPolynomials.coefficients(p)]

    X = gens(R)
    
    Mon = [prod(X[k]^e[k] for k in 1:length(e)) for e in Exp]

    pol = sum(Coeff[i]*Mon[i] for i in 1:length(Coeff))

end

function convert_coeff(p, C::Type=Float64)

    Mon = (monomials(p))
    Coeff = [C(c) for c in coefficients(p)]
    
    pol = sum(Coeff[i]*Mon[i] for i in 1:length(Coeff))

    return pol
end

function convert_DP(P::Vector, C::Type=Float64)
    R = parent(P[1])
    n = length(AbstractAlgebra.gens(R))
    X = (DynamicPolynomials.@polyvar x[1:n] monomial_order = Graded{Reverse{LexOrder}})[1]

    #C = typeof(first(AbstractAlgebra.coefficients(P[1])))

    return [as_polynomial(p, X, C) for p in P]
end

