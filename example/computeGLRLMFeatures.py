import itk, sys

if len(sys.argv) != 9:
    print("Usage: " + sys.argv[0] + " <inputImagePath> <maskImagePath>"
                                    " <NumberOfBinsPerAxis> <PixelValueMin> "
                                    "<PixelValueMax> <DistanceValueMin> "
                                    "<DistanceValueMax> <NeighborhoodRadius>")
    sys.exit(1)

im = itk.imread(sys.argv[1])
mask = itk.imread(sys.argv[2])

filtr = itk.RunLengthTextureFeaturesImageFilter.New(im)
filtr.SetMaskImage(mask)
filtr.SetNumberOfBinsPerAxis(int(sys.argv[3]))
filtr.SetSetHistogramValueMinimum(int(sys.argv[4]))
filtr.SetSetHistogramValueMaximum( int(sys.argv[5]))
filtr.SetHistogramDistanceMinimum(float(sys.argv[6]))
filtr.SetHistogramDistanceMaximum(float(sys.argv[7]))
filtr.SetNeighborhoodRadius([int(sys.argv[8]),int(sys.argv[8]),int(sys.argv[8])])

result = filtr.GetOutput()

itk.imwrite(result, "result.nrrd")
