#include "itkRunLengthTextureFeaturesImageFilter.h"

#include "itkImage.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhood.h"

int main(int argc, char * argv[])
{
    if( argc != 10 )
    {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFilePath> <MaskFilePath> <OutputFilePath> ";
    std::cerr << " <NumberOfBinsPerAxis> <PixelValueMin> ";
    std::cerr << " <PixelValueMax>  <DistanceValueMin> ";
    std::cerr << " <DistanceValueMax> <NeighborhoodRadius> ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

    // Setup types
    using InputImageType = itk::Image< int, 3 >;
    using OutputImageType = itk::Image< itk::Vector< float, 10 > , 3 >;
    using readerType = itk::ImageFileReader< InputImageType >;
    using NeighborhoodType = itk::Neighborhood<typename InputImageType::PixelType,
      InputImageType::ImageDimension>;
    NeighborhoodType neighborhood;

    // Create and setup a reader
    readerType::Pointer reader = readerType::New();
    reader->SetFileName( argv[1] );

    // Create and setup a maskReader
    readerType::Pointer maskReader = readerType::New();
    maskReader->SetFileName( argv[2] );

    // Apply the filter
    using FilterType = itk::Statistics::RunLengthTextureFeaturesImageFilter
                                < InputImageType, OutputImageType >;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());
    filter->SetMaskImage(maskReader->GetOutput());
    filter->SetNumberOfBinsPerAxis(std::atoi(argv[4]));
    filter->SetHistogramValueMinimum(std::atof(argv[5]));
    filter->SetHistogramValueMaximum(std::atof(argv[6]));
    filter->SetHistogramDistanceMinimum(std::atof(argv[7]));
    filter->SetHistogramDistanceMaximum(std::atof(argv[8]));
    neighborhood.SetRadius( std::atoi(argv[9]) );
    filter->SetNeighborhoodRadius(neighborhood.GetRadius());

    // Create and setup a writter
    using WriterType = itk::ImageFileWriter< OutputImageType  >;
    WriterType::Pointer writer = WriterType::New();
    std::string outputFilename = argv[3];
    writer->SetFileName(outputFilename.c_str());
    writer->SetInput(filter->GetOutput());
    writer->Update();

  return EXIT_SUCCESS;
}
