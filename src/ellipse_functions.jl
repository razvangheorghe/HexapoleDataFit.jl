using ImageFiltering
using Images

export rotation_matrix, ellipse_to_points, kernelize, x_bounds_from_ellipse, stdev_weights, get_x_given_y, get_y_given_x

function rotation_matrix(angle::Real; ccw=true)
    #= Helper for rotating points in 2d space

    Builds the rotation matrix specified by the angle and the direction of rotation

    Args :
        angle : angle in radians of rotation
        ccw : whether the angle is measured counter clockwise wrt the positive x axis

    Returns :
        A 2x2 rotation matrix
    =#

    rotate_mat = zeros(2,2);
    if ccw
        rotate_mat = [cos(angle) -sin(angle);
                      sin(angle) cos(angle)];
    else
        rotate_mat = [cos(angle) sin(angle);
                      -sin(angle) cos(angle)];
    end

    return rotate_mat;
end

# might add method for ConicFormEllipse or unite both Conic and Parametric in another structure
function ellipse_to_points(ellipse::ParametricFormEllipse; n=1000::Int)

    theta_plot_vals = range(0, 2*pi, length=n)

    unit_circle = [cos.(theta_plot_vals) sin.(theta_plot_vals)]

    onaxis_ellipse = reshape(ellipse.axis_lengths ./ 2, (1,2)) .* unit_circle

    U = rotation_matrix(deg2rad(ellipse.ccw_angle))

    rotated_ellipse = (U * onaxis_ellipse')'

    shifted_ellipse = vec(ellipse.center)' .+ rotated_ellipse

    return shifted_ellipse
end

function kernelize(ellipse::ParametricFormEllipse) #; gaussian = true)
    X = round.(Int, ellipse_to_points(ellipse, n = 1000))

# TODO change translateImage so that the new image is centered at x(x, y)
    ell = centered(zeros(Float64, (301, 301)))
    for i ∈ 1:size(X, 1)
        ell[X[i, 2], -X[i, 1]] = 1.0
    end
    ell = parent(ell)

    # if gaussian
    gauss3x3 = centered(1/16 .* [1 2 1; 2 4 2; 1 2 1])
    ell = imfilter(ell, gauss3x3)
    # end

    slice1 = mapslices(col -> any(col .!= 0), ell, dims = 1)
    slice2 = mapslices(col -> any(col .!= 0), ell, dims = 2)


    return centered(ell[:, vec(slice1)][vec(slice2), :])
end


function x_bounds_from_ellipse(ellipse::ConicFormEllipse)
    #=
    general conic form ellipse:
    Ax² + Bxy + Cy² + Dx + Ey + F = 0 (1)

    to find extremum points, take ∂/∂y of both sides
    ∴ Bx + 2Cy + E = 0 → get y as f(x) and plug it back into (1)

    → get a quadratic equation in x: ax² + bx + c
    which gives the leftmost and rightmost indices of x as solution
    =#
    A = ellipse.A; B = ellipse.B; C = ellipse.C; D = ellipse.D; E = ellipse.E; F = ellipse.F

    # solving this quadratic gives the x coordinates of the top and bottom
    # they arise when the ∂/∂x is 0
    # a = -A + 4 * A^2 * C / B^2
    # b = 4 * A * C * D / B^2 - 2 * A * E / B
    # c = C * D^2 + B^2 - D * E / B + F

    a = A - B^2 / (4 * C)
    b = D - B * E / (2 * C)
    c = F - E^2 / (4 * C)


    Δ = b^2 - 4 * a * c

    leftmost_x_ind = floor(Int, (-b + √Δ) / (2 * a)) # ceil of these will be points just outside of the ellipse
    rightmost_x_ind = ceil(Int, (-b - √Δ) / (2 * a)) # floor of these will be points just outside of the ellipse

    return leftmost_x_ind, rightmost_x_ind
end

function x_bounds_from_ellipse(ellipse::ParametricFormEllipse)
    ellipse = parametric2conic(ellipse)
    x_bounds_from_ellipse(ellipse)
end



function stdev_weights(ellipse::ConicFormEllipse)

    # find leftmost and rightmost x coordinates of the ellipse
    right_col_ind, left_col_ind = x_bounds_from_ellipse(ellipse)

    xs = Int.([left_col_ind:right_col_ind...])

    # find topmost and bottommost y coordinates of the ellipse
    # for each x between the leftmost and rightmost x coordinate
    ys = [get_y_given_x(col_ind, ellipse) for col_ind ∈ xs]
    ys = StructArray(ys)

    # these are the arrays of topmost and bottommost y coordinates respectively
    top_ys, bottom_ys = fieldarrays(ys) # floats
    top_ys = OffsetArray(ceil.(Int, top_ys), left_col_ind:right_col_ind) # integers
    bottom_ys = OffsetArray(floor.(Int, bottom_ys), left_col_ind:right_col_ind) # integers

    # initialise weights matrix
    weights = centered(zeros(Float64, (301, 301)))

    for x ∈ xs
        # for each column, if the points are inside the ellipse, make them 1.0
        weights[0 - top_ys[x]:0 - bottom_ys[x], x] = ones(length(bottom_ys[x]:top_ys[x]))
    end

    return weights
end

function stdev_weights(ellipse::ParametricFormEllipse)
    ellipse = parametric2conic(ellipse)
    stdev_weights(ellipse)
end

function get_y_given_x(x::Real, ellipse::ConicFormEllipse)
    #=
    general conic form ellipse:
    Ax² + Bxy + Cy² + Dx + Ey + F = 0 (1)

    knowing x, (1) turns into a quadratic in y with two solutions
    ay² + by + c = 0

    where
    a = C
    b = Bx + E
    c = Ax² + Dx + F
    =#
    A = ellipse.A; B = ellipse.B; C = ellipse.C; D = ellipse.D; E = ellipse.E; F = ellipse.F

    a = C
    b = B * x + E
    c = A * x^2 + D * x + F

    Δ = b^2 - 4 * a * c

    y1 = (-b + sqrt(Δ)) / (2 * a)
    y2 = (-b - sqrt(Δ)) / (2 * a)

    return y1, y2
end

function get_x_given_y(y::Real, ellipse::ConicFormEllipse)
    #=
    general conic form ellipse:
    Ax² + Bxy + Cy² + Dx + Ey + F = 0 (1)

    knowing y, (1) turns into a quadratic in x with two solutions
    ax² + xy + c = 0

    where
    a = A
    b = By + D
    c = Cy² + Ey + F
    =#
    A = ellipse.A; B = ellipse.B; C = ellipse.C; D = ellipse.D; E = ellipse.E; F = ellipse.F

    a = C
    b = B * y + D
    c = C * y^2 + E * y + F

    Δ = b^2 - 4 * a * c

    x1 = (-b + sqrt(Δ)) / (2 * a)
    x2 = (-b - sqrt(Δ)) / (2 * a)

    return x1, x2
end
