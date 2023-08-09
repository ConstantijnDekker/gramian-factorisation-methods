import LinearAlgebra as LA

function polar(A)
    (U, _, V) = LA.svd(A)
    U * V'
end

function alterProj(C, tolerance, maxiters)
    Q = LA.I
    X = Q*C
    iters = 0
    #error = max(0, -minimum(X))
    while minimum(X) < -tolerance && iters < maxiters
        Y = map(x -> max(x, 0), X)
        Q = polar(Y*C')
        X = Q*C
        iters += 1
        #push!(error, max(0, -minimum(X)))
    end
    Q, iters#, error
end
