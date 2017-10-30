import itk, sys

if len(sys.argv) != 9:
    print("Usage: " + sys.argv[0] + " <inputImagePath> <maskImagePath>"
                                    " <NumberOfBinsPerAxis> <PixelValueMin> "
                                    "<PixelValueMax> <DistanceValueMin> "
                                    "<DistanceValueMax> <NeighborhoodRadius>")
    sys.exit(1)


Dimension = 3

#Input scan reader
InputPixelType = itk.ctype('signed short')
InputImageType = itk.Image[InputPixelType, Dimension]
imReader = itk.ImageFileReader[InputImageType].New()
imReader.SetFileName(sys.argv[1])

#Input mask reader
MaskPixelType = itk.ctype('unsigned char')
MaskImageType = itk.Image[MaskPixelType, Dimension]
maskReader = itk.ImageFileReader[MaskImageType].New()
maskReader.SetFileName(sys.argv[2])

im = imReader.GetOutput()
mask = maskReader.GetOutput()

filtr = itk.RunLengthTextureFeaturesImageFilter.New(im)
filtr.SetMaskImage(mask)
filtr.SetNumberOfBinsPerAxis(int(sys.argv[3]))
filtr.SetHistogramValueMinimum(int(sys.argv[4]))
filtr.SetHistogramValueMaximum( int(sys.argv[5]))
filtr.SetHistogramDistanceMinimum(float(sys.argv[6]))
filtr.SetHistogramDistanceMaximum(float(sys.argv[7]))
filtr.SetNeighborhoodRadius([int(sys.argv[8]),int(sys.argv[8]),int(sys.argv[8])])

result = filtr.GetOutput()

itk.imwrite(result, "result.nrrd")
