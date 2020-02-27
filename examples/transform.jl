using HexapoleDataFit
using CoordinateTransformations

img = readImage("averaged_V+H_NO.dat")
rotimg = rotateImage(img, -102)
trnsimg = translateImage(rotimg, (50, 0))
trnssimg = translateImage(trnsimg, (0, 50))

showImage(translateImage(rotimg, CartesianIndex(0, 50)))
showImage(translateImage(rotimg, CartesianIndex(50, 0)))



# show (in Juno) rotated, original and translated image
showImage(img)
showImage(rotimg)
showImage(trnsimg)
showImage(trnssimg)


# add circles over an image
showImage(img)
plotRadius(90)
plotRadius(112)
plotRadius(130)

# add ellipses over an image
showImage(trnsimg)
ellipse = ParametricFormEllipse([240.0 220.0], [0.0 0.0], 131)
plotEllipse(ellipse)
