using LinearAlgebra

export fitEllipse

function fitEllipse(x::Array{Float64, 1}, y::Array{Float64, 1})
    # outputs [A; B; C; D; E; F], a Vector containing the parameters of the equation
    # Ax² + Bxy + Cy² + Dx + Ey + F = 0
    #
    # Based on Direct Least Squares Fitting of Ellipses
    # by A. W. Fitzgibbon, M. Pilu, R. B. Fisher
    # http://cseweb.ucsd.edu/~mdailey/Face-Coord/ellipse-specific-fitting.pdf


    # Build design matrix
    D = [x.*x x.*y y.*y x y ones(Float64, size(x))]

    # Build scatter matrix
    S = transpose(D) * D

    # Build 6x6 constraint matrix
    C = zeros(Float64, (6, 6))
    C[1, 3] = C[3, 1] = -2
    C[2, 2] = 1

    # Solve generalised eigensystem
    geval, gevec = eigen(S, C)
    # geval - eigenvalues
    # gevec - eigenvectors

    # Find the only negative eigenvalue index
    negindex = findall(x->x<0, geval)

    # Get fitted parameters
    return a = gevec[:, negindex]

end
