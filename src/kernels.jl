using ImageFiltering
using Images

export rotation_matrix, ellipse_to_points, kernelize

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


# i'll probably never use this
function kernelize(ellipse_image::OffsetArray{Float64,2,Array{Float64,2}})
    ell = parent(ellipse_image)
    slice1 = mapslices(col -> any(col .!= 0), ell, dims = 1)
    slice2 = mapslices(col -> any(col .!= 0), ell, dims = 2)

    return centered(ell[:, vec(slice1)][vec(slice2), :])
end
