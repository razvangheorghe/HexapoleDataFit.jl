# adapted from https://github.com/chaaland/EllipseFit.jl by Casey Haaland

export ConicFormEllipse, ParametricFormEllipse

abstract type EllipseForm end

struct ConicFormEllipse{T<:Real} <: EllipseForm
    A::T
    B::T
    C::T
    D::T
    E::T
    F::T

    function ConicFormEllipse(A::T, B::T, C::T, D::T, E::T, F::T) where {T<:Real}
        if B ^2 - 4 * A * C >= 0
            error("Discriminant is non-negative. Input does not denote an ellipse")
        end
    end
end
ConicFormEllipse(A::Real, B::Real, C::Real, D::Real, E::Real, F::Real) = ConicFormEllipse(promote(A, B, C, D, E, F)...)


struct ParametricFormEllipse{T<:Real, U<:Real, V<:Real}
    axis_lengths::Array{T} # first element is width, second element is height
    center::Array{U} # first element is x₀ coordinate, second element is y₀ coordinate
    ccw_angle::V # counterclockwise angle of rotation of the ellipse in degrees

    function ParametricFormEllipse(axis_lengths::Array{T}, center::Array{U}, ccw_angle::V) where {T<:Real, U<:Real, V<:Real}
        axis_lengths = vec(axis_lengths)
        center = vec(center)

        if size(center)[1] != 2
            error("Parameter 'center' must be vector of length 2")
        elseif size(axis_lengths)[1] != 2
            error("Parameter 'axis_lengths' must be vector of length 2")
        elseif sum(axis_lengths .< 0) > 0
            error("Parameter 'axis_lengths' must be nonnegative")
        end
    end
end



function conic2parametric(cform::ConicFormEllipse)
    # TODO reference the different equations properly

    # mostly from here: https://en.wikipedia.org/wiki/Ellipse#General_ellipse

    #=
    Helper for converting from conic form to parametric form of ellipse
    Given an ellipse in standard conic section form
            A * x^2 + B * x * y + C * y^2 + D * x + E * y + F = 0
    convert it to parametric form
            [x_c y_c] + rot_mat2d(ccw_angle) * [a*cos(theta) b*sin(theta)]
    Args :

    Returns :
        semiaxis_lengths : array of the form [a b] where a is half the length
                           of the axis aligned ellipse along the x-axis and b
                           is half the length along the y-axis (before rotation)
        center : array of the x and y coordinates of the center of the ellipse
        ccw_angle : The counter clockwise angle (in degrees) to rotate the ellipse wrt
                    the positive x-axis
    =#

    A = cform.A; B = cform.B; C = cform.C; D = cform.D; E = cform.E; F = cform.F

    A_Q = [A B/2 D/2; B/2 C E/2; D/2 E/2 F]
    Δ = det(A_Q) # -1/4 * (AE² + CD² + FB² -BDE - 4 * ACF)

    A_33 = [A B/2; B/2 C]
    J = det(A_33) # -1/4 * (B² - 4 * AC)

    I = tr(A_33)

    sq = sqrt((A - C)^2 + B^2)

    #=
    now write the equations (21) and (22) at http://mathworld.wolfram.com/Ellipse.html
    in terms of Δ, J, I and sq --> semiaxes a and b
    =#

    a = sqrt(-2 * Δ * (I + sq)) / (2*J)
    b = sqrt(-2 * Δ * (I - sq)) / (2*J)

    # a = 2 * sqrt(2 * Δ / (J * (sq - I)))
    # b = 2 * sqrt(2 * Δ / (J * (-sq - I)))

    # a = -2 * sqrt(2 * (A*E^2 + C*D^2 - B*D*E + (B^2 - 4*A*C)*F) * ((A + C) + (sqrt((A - C)^2 + B^2)))) / (B^2 - 4*A*C)
    # b = -2 * sqrt(2 * (A*E^2 + C*D^2 - B*D*E + (B^2 - 4*A*C)*F) * ((A + C) - (sqrt((A - C)^2 + B^2)))) / (B^2 - 4*A*C)

    axis_lengths = [2*a 2*b]

    x₀, y₀ = inv(A_33) * [-D/2; -E/2]

    # x₀ = (2 * C * D - B * E) / ( -16 * J)
    # y₀ = (2 * A * E - B * D) / ( -16 * J)

    center = [x₀ y₀]

    # θ = atan(1/B * (C - A - sq))

    # θ is in radians
    if A <= C
        θ = 1/2 * acot((A - C) / (2 * B))
    else
        θ = π / 2 + 1/2 * acot((A - C) / (2 * B))
    end

    ccw_angle = rad2deg(θ)

    return ParametricFormEllipse(axis_lengths, center, ccw_angle)
end
