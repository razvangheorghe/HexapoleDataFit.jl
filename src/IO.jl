using DelimitedFiles
using OffsetArrays
using PyPlot
using ImageFiltering


export readImage, readSimImage, writeImage, saveImage

function readImage(path::String)
        image = readdlm(path, Float64)
        z = reshape(image[:, 3], (301, 301))
        z = centered(z) # creates an OffsetArray where the x and y indices are from -150 to 150 instead of 1 to 301
        return z
end

function readSimImage(path::String, column_nr::Int64)
        image = readdlm(path, Float64)
        z = reshape(image[:, column_nr], (301, 301))
        z = copy(transpose(centered(z)))
        return z
end

function writeImage(path::String, image::OffsetArray{Float64,2,Array{Float64,2}})
        x = reshape(repeat(transpose(collect(axes(image)[1])), 301, 1), (90601, 1))
        y = reshape(repeat(collect(reverse(axes(image)[2])), 1, 301), (90601, 1))
        z = reshape(image, (90601, 1))
        open(path, "w") do io
           writedlm(io, [x y z])
        end
end

function saveImage(path::String, image::Array{Float64,2})
# save current figure as pdf (path = name.pdf) or as jpg (path = name.jpg )
        savefig(path, bbox_inches="tight")
end

function saveImage(path::String, image::OffsetArray{Float64,2,Array{Float64,2}})
# save current figure as pdf (path = name.pdf) or as jpg (path = name.jpg )
        savefig(path, bbox_inches="tight")
end
