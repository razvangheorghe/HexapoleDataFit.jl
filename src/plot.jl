using PyPlot
using OffsetArrays

export showImage, plotRadius, plotEllipse

function showImage(image::OffsetArray{Float64,2,Array{Float64,2}})
        clf()
        imshow(image, aspect = "equal", extent = (minimum(axes(image)[1]), maximum(axes(image)[1]), minimum(axes(image)[2]), maximum(axes(image)[2])))
        colorbar()
        gcf()
end

function showImage(image::Array{Float64,2})
        clf()
        imshow(image, aspect = "equal", extent = (minimum(axes(image)[1]), maximum(axes(image)[1]), minimum(axes(image)[2]), maximum(axes(image)[2])))
        colorbar()
        gcf()
end

function showImage(x::Array{Float64,2}, y::Array{Float64,2}, z::Array{Float64,2})
        clf()
        pcolormesh(x, y, z)
        xlim([minimum(x), maximum(x)]);
        ylim([minimum(y), maximum(y)]);
        gca().set_aspect("equal")
        colorbar()
        gcf()
end


function plotRadius(radius::Prec64)
        plt.gca().add_artist(plt.Circle((0.0, 0.0), radius = radius, color = "r", fill=false))
        gcf()
end

function plotEllipse(center::Tuple{Float64,Float64}, width::Float64, height::Float64, angle::Float64)
        ax = gca()
        ell = matplotlib.patches.Ellipse(center, width, height, angle, color = "r", fill = false)
        ax.add_patch(ell)
        gcf()
end


function plotEllipse(ell::ParametricFormEllipse)
        ax = gca()
        center = tuple(ell.center...)
        width = ell.axis_lengths[1]
        height = ell.axis_lengths[2]
        angle = ell.ccw_angle
        ell = matplotlib.patches.Ellipse(center, width, height, angle, color = "r", fill = false)
        ax.add_patch(ell)
        gcf()
end


function plotEllipse(ell::ConicFormEllipse)
        plotEllipse(conic2parametric(ell))
end
