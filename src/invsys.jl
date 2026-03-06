export invsys

function invsys(s::Series{C, M}, X = variables(s), d::Int = maxdegree(s)) where {C, M <: AbstractMonomial }

    D = monomials(X,0:d)
    Ip = Series{C,M}[]
    z0 = zero(s)
    for m in D
        sn = m*s
        if sn != z0
            push!(Ip,sn)
        end
    end
    Ip
end

