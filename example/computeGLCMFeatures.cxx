#include "itkCoocurrenceTextureFeaturesImageFilter.h"

#include "itkImage.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhood.h"

int
main(int argc, char * argv[])
{
  if (argc != 8)
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFilePath> <MaskFilePath> <OutputFilePath> ";
    std::cerr << " <NumberOfBinsPerAxis> <PixelValueMin> ";
    std::cerr << " <PixelValueMax> <NeighborhoodRadius> ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  // Setup types
  using InputImageType = itk::Image<int, 3>;
  using OutputImageType = itk::Image<itk::Vector<float, 10>, 3>;
  using readerType = itk::ImageFileReader<InputImageType>;
  using NeighborhoodType = itk::Neighborhood<typename InputImageType::PixelType, InputImageType::ImageDimension>;
  NeighborhoodType neighborhood;

  // Create and setup a reader
  readerType::Pointer reader = readerType::New();
  reader->SetFileName(argv[1]);

  // Create and setup a maskReader
  readerType::Pointer maskReader = readerType::New();
  maskReader->SetFileName(argv[2]);

  // Apply the filter
  using FilterType = itk::Statistics::CoocurrenceTextureFeaturesImageFilter<InputImageType, OutputImageType>;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(reader->GetOutput());
  filter->SetMaskImage(maskReader->GetOutput());
  filter->SetNumberOfBinsPerAxis(std::stoi(argv[4]));
  filter->SetHistogramMinimum(std::stod(argv[5]));
  filter->SetHistogramMaximum(std::stod(argv[6]));
  neighborhood.SetRadius(std::stoi(argv[7]));
  filter->SetNeighborhoodRadius(neighborhood.GetRadius());

  // Create and setup a writter
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  std::string         outputFilename = argv[3];
  writer->SetFileName(outputFilename.c_str());
  writer->SetInput(filter->GetOutput());
  writer->Update();

  return EXIT_SUCCESS;
}
