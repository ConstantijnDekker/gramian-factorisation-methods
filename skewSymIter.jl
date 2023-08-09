import LinearAlgebra as LA
import JuMP
import Tulip

# Return skewsymmetric S such that (I + S) Y >= 0
function g(Y)
    find_skew(Y, Tulip)
end

# Find skew-symmetric solution to (I + S) Y >= 0 using OPT module.
function find_skew(Y, Opt)
    (m, _) = size(Y)
    model = JuMP.Model(Opt.Optimizer)
    # Make sure the model does not output unnecessary status information.
    JuMP.set_silent(model)

    JuMP.@variable(model, S[1:m, 1:m] in JuMP.SkewSymmetricMatrixSpace())

    JuMP.@constraint(model, (LA.I + S) * Y .>= 0)

    # t models the maximum absolute value of S
    JuMP.@variable(model, t)
    JuMP.@constraint(model, t .>= S)
    JuMP.@objective(model, Min, t)

    JuMP.optimize!(model)

    # Return obtained solutions and whether it was succesful
    return (JuMP.value.(S), JuMP.termination_status(model) == JuMP.OPTIMAL)
end

# Select orthogonal approximation of I + S
function f(S, k)
    if k == 3
        (LA.I + S) * LA.inv(LA.sqrt(LA.I - S^2))
    elseif k == 2
        LA.exp(S)
    elseif k == 1
        S + LA.sqrt(LA.I + S^2)
    end
end

function skewSymIter(C, k, tolerance, maxiters)
    (m, n) = size(C)
    Q = Matrix{eltype(C)}(LA.I, m, m)
    Y = Q * C
    iters = 0
    while minimum(Y) < -tolerance && iters < maxiters
        S, feasible = g(Y)
        V = if !feasible || LA.opnorm(S) > 1
            Z = map(x -> max(x, 0), Y)
            polar(Z * Y')
        else
            f(S, k)
        end
        Y = V * Y
        Q = V * Q
        iters += 1
    end

    return (Q, iters)
end
