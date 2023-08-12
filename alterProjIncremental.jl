import LinearAlgebra as LA

function polar(A)
    (U, _, V) = LA.svd(A)
    U * V'
end

function alterProjIncremental(C, tolerance, maxiters)
    Q = LA.I
    V = LA.I
    X = Q*C
    iters = 0
    #error = [max(0, -minimum(X))]
    while minimum(X) < -tolerance && iters < maxiters
        Y = map(x -> max(x, 0), X)
        V = polar(Y * X')
        X = V*X
        Q = V*Q
        iters += 1
        #push!(error, max(0, -minimum(X)))
    end
    Q, iters#, error
end
