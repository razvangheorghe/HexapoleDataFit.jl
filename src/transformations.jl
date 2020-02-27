using ImageCore
using ImageTransformations
using CoordinateTransformations
using OffsetArrays
using Images

export translateImage, rotateImage, warp2translate_cartesian, ellipse2circle_image

function warp2translate_cartesian(cart::CartesianIndex{2})
    return CartesianIndex(-cart[2], cart[1])
end

function translateImage(image::OffsetArray{Float64,2,Array{Float64,2}}, xyshift::Tuple{Prec64,Prec64})
    xyshift = tuple(xyshift[2], -xyshift[1])

    t = Translation(xyshift)
    translated_image = warp(image, t, indices_spatial(image), 0.0)

    return translated_image
end

function translateImage(image::OffsetArray{Float64,2,Array{Float64,2}}, xyshift::CartesianIndex{2})
    x₀ = xyshift[2]
    y₀ = -xyshift[1]
    t = Translation((x₀, y₀))
    translated_image = warp(image, t, indices_spatial(image), 0.0)

    return translated_image
end

function translateImage(image::Array{Float64,2}, xyshift::Tuple{Prec64,Prec64})
    xyshift = tuple(xyshift[2], -xyshift[1])

    t = Translation(xyshift)
    translated_image = warp(image, t, indices_spatial(image), 0.0)

    return translated_image
end

function rotateImage(image::OffsetArray{Float64,2,Array{Float64,2}}, θ::Prec64)
    θ = deg2rad(θ)
    rotated_image = imrotate(image, -θ, indices_spatial(image), 0.0)

    return rotated_image
end

function rotateImage(image::Array{Float64,2}, θ::Prec64)
    θ = deg2rad(θ)
    rotated_image = imrotate(image, -θ, indices_spatial(image), 0.0)

    return rotated_image
end


function ellipse2circle_image(image::OffsetArray{Float64,2,Array{Float64,2}}, ellipse::ParametricFormEllipse)
# takes in a centered stretched ion image and outputs a squished, centered circular ion image
# and the radius of the best fit circle
# currently the radius of the new image is the average of the ellipse height and width

    width = ellipse.axis_lengths[1]
    height = ellipse.axis_lengths[2]
    ellipse_angle = ellipse.ccw_angle

# rotate image such that the ellipse's axes are allilgned with x and y
    image2 = rotateImage(image, -ellipse_angle)

# define the radius of the final circle - average of width and height in this case
    diameter = (width + height) / 2
    x_scale = diameter / width
    y_scale = diameter / height

# scale x axis
    image3 = centered(imresize(image2, round(Int, 301 * x_scale), size(image2)[2]))

# scale y axis
    image4 = centered(imresize(image3, size(image3)[1], round(Int, 301 * y_scale)))

# fit inside new 301 x 301 frame
    image5 = centered(zeros(Float64, (301, 301)))
    if length(indices_spatial(image4)[1]) > 301
        if length(indices_spatial(image4)[2]) > 301
            image5[:, :] = image4[-150:150, -150:150]
        else
            image5[:, indices_spatial(image4)[2]] = image4[-150:150, :]
        end
    else
        if length(indices_spatial(image4)[2]) > 301
            image5[indices_spatial(image4)[1], :] = image4[:, indices_spatial(image5)[2]]
        else
            image5[indices_spatial(image4)[1], indices_spatial(image4)[2]] = image4[:, :]
        end
    end
# the result here is kind of weird. TODO check again img5 in test_files\\flatten_ellipse.jl

# need to recenter again due to pixels being discrete
    circle = kernelize(ParametricFormEllipse([width*x_scale height*y_scale], [0.0 0.0], 0.0))
    maxvalue, maxcoord = findmax(imfilter(image5, circle))
    maxcoord = warp2translate_cartesian(maxcoord)
# actual recentering
    image6 = translateImage(image5, maxcoord)

# rotate back
    return rotateImage(image6, ellipse_angle), diameter/2
end
