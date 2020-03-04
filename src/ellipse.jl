# adapted from https://github.com/chaaland/EllipseFit.jl by Casey Haaland

using LinearAlgebra

export ParametricFormEllipse, ConicFormEllipse, parametric2conic, conic2parametric



abstract type EllipseForm end

struct ParametricFormEllipse{T<:Real, U<:Real, V<:Real}
    axis_lengths::Array{T}
    center::Array{U}
    ccw_angle::V

    function ParametricFormEllipse(axis_lengths::Array{T}, center::Array{U}, ccw_angle::V) where {T<:Real, U<:Real, V<:Real}
        axis_lengths = vec(axis_lengths)
        center = vec(center)

        if size(center)[1] != 2
            error("Parameter 'center' must be vector of length 2")
        elseif size(axis_lengths)[1] != 2
            error("Parameter 'axis_lengths' must be vector of length 2")
        elseif sum(axis_lengths .< 0) > 0
            error("Parameter 'axis_lengths' must be nonnegative")
        else
            new{T,U,V}(axis_lengths, center, ccw_angle)
        end
    end
end

function parametric2conic(ellipse::ParametricFormEllipse)
    # inputs
    a = ellipse.axis_lengths[1] / 2
    b = ellipse.axis_lengths[2] / 2
    x₀ = ellipse.center[1]
    y₀ = ellipse.center[2]
    θ = deg2rad(ellipse.ccw_angle)

    # outputs
    A = a^2 * sin(θ)^2 + b^2 * cos(θ)^2
    B = 2 * (b^2 - a^2) * sin(θ) * cos(θ)
    C = a^2 * cos(θ)^2 + b^2 * sin(θ)^2
    D = -2 * A * x₀ - B * y₀
    E = -B * x₀ - 2 * C * y₀
    F = A * x₀^2 + B * x₀ * y₀ + C * y₀^2 - a^2 * b^2

    return ConicFormEllipse(A, B, C, D, E, F)
end

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
        new{T}(A, B, C, D, E, F)
    end
end

ConicFormEllipse(A::Real, B::Real, C::Real, D::Real, E::Real, F::Real) = ConicFormEllipse(promote(A, B, C, D, E, F)...)


function conic2parametric(ellipse::ConicFormEllipse)
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
        axis_lengths : array of the form [2*a 2*b] where a is half the width
                          and b is half the height (before rotation)
        center : array of the x and y coordinates of the center of the ellipse
        ccw_angle : The counter clockwise angle (in degrees) to rotate the ellipse wrt
                    the positive x-axis
    =#

    A = ellipse.A; B = ellipse.B; C = ellipse.C; D = ellipse.D; E = ellipse.E; F = ellipse.F

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

    axis_lengths = [2*a 2*b]

    x₀, y₀ = inv(A_33) * [-D/2; -E/2]

    center = [x₀ y₀]

    # all of the angles are mapped onto [-90,90]
    # as the system has a period of 180°
    # i.e. [-280, -100, 80, 260, etc] → 80

    if B == 0
        if A <= C
            θ = 0
        else
            θ = π/2
        end
    else
        θ = atan(1/B * (C - A - sq))
    end

    ccw_angle = rad2deg(θ)


    return ParametricFormEllipse(axis_lengths, center, ccw_angle)
end
