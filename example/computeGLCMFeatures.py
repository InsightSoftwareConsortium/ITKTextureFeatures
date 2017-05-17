import itk, sys

if len(sys.argv) != 7:
    print("Usage: " + sys.argv[0] + " <inputImagePath> <maskImagePath>"
                                    " <NumberOfBinsPerAxis> <PixelValueMin> "
                                    "<PixelValueMax> <NeighborhoodRadius>")
    sys.exit(1)

im = itk.imread(sys.argv[1])
mask = itk.imread(sys.argv[2])

filtr = itk.ScalarImageToTextureFeaturesImageFilter.New(im)
filtr.SetMaskImage(mask)
filtr.SetNumberOfBinsPerAxis(int(sys.argv[3]))
filtr.SetPixelValueMinMax(int(sys.argv[4]), int(sys.argv[5]))
filtr.SetNeighborhoodRadius([int(sys.argv[6]),int(sys.argv[6]),int(sys.argv[6])])

result = itk.VectorImage[itk.F,3].New()
result = filtr.GetOutput()

itk.imwrite(result, "result.nrrd")
