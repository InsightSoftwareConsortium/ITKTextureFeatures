/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkCoocurrenceTextureFeaturesImageFilter.h"

#include "itkImage.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhood.h"
#include "itkTestingMacros.h"

int
CoocurrenceTextureFeaturesImageFilterTestWithoutMask(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Missing parameters." << std::endl;
    std::cerr << "Usage: " << argv[0] << " inputImageFile"
              << " outputImageFile"
              << " [numberOfBinsPerAxis]"
              << " [pixelValueMin]"
              << " [pixelValueMax]"
              << " [neighborhoodRadius]" << std::endl;
    return EXIT_FAILURE;
  }

  constexpr unsigned int ImageDimension = 3;
  constexpr unsigned int VectorComponentDimension = 8;

  // Declare types
  using InputPixelType = float;
  using OutputPixelComponentType = float;
  using OutputPixelType = itk::Vector<OutputPixelComponentType, VectorComponentDimension>;

  using InputImageType = itk::Image<InputPixelType, ImageDimension>;
  using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;
  using ReaderType = itk::ImageFileReader<InputImageType>;
  using NeighborhoodType = itk::Neighborhood<InputImageType::PixelType, InputImageType::ImageDimension>;

  // Create and set up a reader
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);

  // Create the filter
  using FilterType = itk::Statistics::CoocurrenceTextureFeaturesImageFilter<InputImageType, OutputImageType>;
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput(reader->GetOutput());

  if (argc >= 4)
  {
    unsigned int numberOfBinsPerAxis = std::stoi(argv[3]);
    filter->SetNumberOfBinsPerAxis(numberOfBinsPerAxis);

    FilterType::PixelType pixelValueMin = std::stod(argv[4]);
    FilterType::PixelType pixelValueMax = std::stod(argv[5]);
    filter->SetHistogramMinimum(pixelValueMin);
    filter->SetHistogramMaximum(pixelValueMax);

    NeighborhoodType::SizeValueType neighborhoodRadius = std::stoi(argv[6]);
    NeighborhoodType                hood;
    hood.SetRadius(neighborhoodRadius);
    filter->SetNeighborhoodRadius(hood.GetRadius());
  }

  ITK_TRY_EXPECT_NO_EXCEPTION(filter->Update());

  // Create and set up a writer
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(filter->GetOutput());

  ITK_TRY_EXPECT_NO_EXCEPTION(writer->Update());


  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}
