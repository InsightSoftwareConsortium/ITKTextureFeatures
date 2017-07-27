/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkCoocurrenceTextureFeaturesImageFilter_hxx
#define itkCoocurrenceTextureFeaturesImageFilter_hxx

#include "itkCoocurrenceTextureFeaturesImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkBinaryFunctorImageFilter.h"

namespace itk
{
namespace Statistics
{
template< typename TInputImage, typename TOutputImage>
CoocurrenceTextureFeaturesImageFilter< TInputImage, TOutputImage >
::CoocurrenceTextureFeaturesImageFilter() :
    m_NumberOfBinsPerAxis( itkGetStaticConstMacro( DefaultBinsPerAxis ) ),
    m_Min( NumericTraits<PixelType>::NonpositiveMin() ),
    m_Max( NumericTraits<PixelType>::max() ),
    m_InsidePixelValue( NumericTraits<PixelType>::OneValue() )
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );

  // Mark the "MaskImage" as an optional named input. First it has to
  // be added to the list of named inputs then removed from the
  // required list.
  Self::AddRequiredInputName("MaskImage");
  Self::RemoveRequiredInputName("MaskImage");

  // Set the offset directions to their defaults: half of all the possible
  // directions 1 pixel away. (The other half is included by symmetry.)
  // We use a neighborhood iterator to calculate the appropriate offsets.
  typedef Neighborhood<typename InputImageType::PixelType,
    InputImageType::ImageDimension> NeighborhoodType;
  NeighborhoodType hood;
  hood.SetRadius( 1 );

  // Select all "previous" neighbors that are face+edge+vertex
  // connected to the iterated pixel. Do not include the currentInNeighborhood pixel.
  unsigned int centerIndex = hood.GetCenterNeighborhoodIndex();
  OffsetVectorPointer offsets = OffsetVector::New();
  for( unsigned int d = 0; d < centerIndex; d++ )
    {
    OffsetType offset = hood.GetOffset( d );
    offsets->push_back( offset );
    }
  this->SetOffsets( offsets );
  NeighborhoodType nhood;
  nhood.SetRadius( 2 );
  this->m_NeighborhoodRadius = nhood.GetRadius( );

  this->m_Normalize = false;
}

template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::SetOffset( const OffsetType offset )
{
  OffsetVectorPointer offsetVector = OffsetVector::New();
  offsetVector->push_back( offset );
  this->SetOffsets( offsetVector );
}

template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{

  typename TInputImage::Pointer input = InputImageType::New();
  input->Graft(const_cast<TInputImage *>(this->GetInput()));

  typedef PreProcessingFunctor PPFType;
  PPFType ppf(m_NumberOfBinsPerAxis, m_InsidePixelValue, m_Min, m_Max);

  typedef BinaryFunctorImageFilter< MaskImageType, InputImageType, InputImageType, PPFType> BinaryFunctorType;
  typename BinaryFunctorType::Pointer functorF = BinaryFunctorType::New();
  if (this->GetMaskImage() != ITK_NULLPTR)
    {
    typename TInputImage::Pointer mask = MaskImageType::New();
    mask->Graft(const_cast<TInputImage *>(this->GetMaskImage()));
    functorF->SetInput1(mask);
    }
  else
    {
    functorF->SetConstant1(m_InsidePixelValue);
    }
  functorF->SetInput2(input);
  functorF->SetFunctor(ppf);
  functorF->SetNumberOfThreads(this->GetNumberOfThreads());

  functorF->Update();
  m_DigitalizedInputImage = functorF->GetOutput();
}

template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  // Free internal image
  this->m_DigitalizedInputImage = ITK_NULLPTR;
}


template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputRegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  // Recuperation of the different inputs/outputs
  OutputImageType* outputPtr = this->GetOutput();

  ProgressReporter progress( this,
                             threadId,
                             outputRegionForThread.GetNumberOfPixels() );

  // Creation of the output pixel type
  typename TOutputImage::PixelType outputPixel;

  // Separation of the non-boundary region that will be processed in a different way
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage > boundaryFacesCalculator;
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType >::FaceListType
  faceList = boundaryFacesCalculator( this->m_DigitalizedInputImage, outputRegionForThread, m_NeighborhoodRadius );
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType >::FaceListType::iterator fit = faceList.begin();

  // Declaration of the variables useful to iterate over the all image region
  bool isInImage;
  IndexType firstIndex;
  for ( unsigned int i = 0; i < this->m_NeighborhoodRadius.Dimension; i++ )
    {
    firstIndex[i] = 0;
    }
  outputPixel = outputPtr->GetPixel(firstIndex);
  typename OffsetVector::ConstIterator offsets;

  // Declaration of the variables useful to iterate over the all the offsets
  OffsetType offset;
  unsigned int totalNumberOfFreq;


  vnl_matrix<unsigned int> hist(m_NumberOfBinsPerAxis, m_NumberOfBinsPerAxis);

  // Declaration of the variables useful to iterate over the all neighborhood region
  PixelType currentInNeighborhoodPixelIntensity;

  // Declaration of the variables useful to iterate over the run
  PixelType pixelIntensity( NumericTraits<PixelType>::ZeroValue() );
  OffsetType tempOffset;

  /// ***** Non-boundary Region *****
  for (; fit != faceList.end(); ++fit )
    {
    NeighborhoodIteratorType inputNIt(m_NeighborhoodRadius, this->m_DigitalizedInputImage, *fit );
    typedef itk::ImageRegionIterator< OutputImageType> IteratorType;
    IteratorType outputIt( outputPtr, *fit );

    // Iteration over the all image region
    while( !inputNIt.IsAtEnd() )
      {
      // If the voxel is outside of the mask, don't treat it
      if( inputNIt.GetCenterPixel() < ( - 5) ) //the pixel is outside of the mask
        {
        progress.CompletedPixel();
        ++inputNIt;
        ++outputIt;
        continue;
        }
      // Initialisation of the histogram
      hist.fill(0);

      totalNumberOfFreq = 0;
      // Iteration over all the offsets
      for( offsets = m_Offsets->Begin(); offsets != m_Offsets->End(); ++offsets )
        {
        offset = offsets.Value();
        // Iteration over the all neighborhood region
        for(NeighborIndexType nb = 0; nb<inputNIt.Size(); ++nb)
          {
          // Test if the current voxel is in the mask and is the range of the image intensity specified
          currentInNeighborhoodPixelIntensity =  inputNIt.GetPixel(nb);
          if( currentInNeighborhoodPixelIntensity < 0 )
            {
            continue;
            }

          // Test if the current offset is still pointing to a voxel inside th neighborhood
          tempOffset = inputNIt.GetOffset(nb) + offset;
          if(!(this->IsInsideNeighborhood(tempOffset)))
          {
            continue;
          }

          // Test if the part of the neighborhood pointed by the offset is still part of the image
          if(fit == faceList.begin())
            {
            inputNIt.GetPixel(tempOffset, isInImage);
            if(!isInImage)
              {
              break;
              }
            }

          // Test if the pointed voxel is in the mask and is the range of the image intensity specified
          pixelIntensity = inputNIt.GetPixel(tempOffset);
          if(pixelIntensity< 0 )
            {
            continue;
            }

          // Increase the corresponding bin in the histogram
          totalNumberOfFreq++;
          hist[currentInNeighborhoodPixelIntensity][pixelIntensity]++;
          }
        }
      // Compute the run length features
      this->ComputeFeatures( hist, totalNumberOfFreq, outputPixel);
      outputIt.Set(outputPixel);

      progress.CompletedPixel();
      ++inputNIt;
      ++outputIt;
      }
    }
}

template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
{
  // Call superclass's version
  Superclass::GenerateOutputInformation();

  OutputImageType* output = this->GetOutput();
  // If the output image type is a VectorImage the number of
  // components will be properly sized if before allocation, if the
  // output is a fixed width vector and the wrong number of
  // components, then an exception will be thrown.
  if ( output->GetNumberOfComponentsPerPixel() != 8 )
    {
    output->SetNumberOfComponentsPerPixel( 8 );
    }
}


template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::SetPixelValueMinMax( PixelType min, PixelType max )
{
  if( this->m_Min != min || this->m_Max != max )
    {
    this->m_Min = min;
    this->m_Max = max;
    this->Modified();
    }
}

template<typename TInputImage, typename TOutputImage>
bool
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::IsInsideNeighborhood(const OffsetType &iteratedOffset)
{
  bool insideNeighborhood = true;
  for ( unsigned int i = 0; i < this->m_NeighborhoodRadius.Dimension; ++i )
    {
    int boundDistance = m_NeighborhoodRadius[i] - Math::abs(iteratedOffset[i]);
    if(boundDistance < 0)
      {
      insideNeighborhood = false;
      break;
      }
    }
  return insideNeighborhood;
}

template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::ComputeFeatures( const vnl_matrix<unsigned int> &hist, const unsigned int totalNumberOfFreq,
                   typename TOutputImage::PixelType &outputPixel)
{
    // Now get the various means and variances. This is takes two passes
    // through the histogram.
    double pixelMean;
    double marginalMean;
    double marginalDevSquared;
    double pixelVariance;

    this->ComputeMeansAndVariances(hist,
                                   totalNumberOfFreq,
                                   pixelMean,
                                   marginalMean,
                                   marginalDevSquared,
                                   pixelVariance);

    // Finally compute the texture features. Another pass.
    MeasurementType energy      = NumericTraits< MeasurementType >::ZeroValue();
    MeasurementType entropy     = NumericTraits< MeasurementType >::ZeroValue();
    MeasurementType correlation = NumericTraits< MeasurementType >::ZeroValue();

    MeasurementType inverseDifferenceMoment  = NumericTraits< MeasurementType >::ZeroValue();

    MeasurementType inertia             = NumericTraits< MeasurementType >::ZeroValue();
    MeasurementType clusterShade        = NumericTraits< MeasurementType >::ZeroValue();
    MeasurementType clusterProminence   = NumericTraits< MeasurementType >::ZeroValue();
    MeasurementType haralickCorrelation = NumericTraits< MeasurementType >::ZeroValue();

    double pixelVarianceSquared = pixelVariance * pixelVariance;
    // Variance is only used in correlation. If variance is 0, then
    //   (index[0] - pixelMean) * (index[1] - pixelMean)
    // should be zero as well. In this case, set the variance to 1. in
    // order to avoid NaN correlation.
    if( Math::FloatAlmostEqual( pixelVarianceSquared, 0.0, 4, 2*NumericTraits<double>::epsilon() ) )
      {
      pixelVarianceSquared = 1.;
      }
    const double log2 = std::log(2.0);

    for(unsigned int a = 0; a < m_NumberOfBinsPerAxis; a++)
      {
      for(unsigned int b = 0; b < m_NumberOfBinsPerAxis; b++)
        {
        float frequency = hist[a][b] / (float)totalNumberOfFreq;
        if ( Math::AlmostEquals( frequency, NumericTraits< float >::ZeroValue() ) )
          {
          continue; // no use doing these calculations if we're just multiplying by
                    // zero.
          }

      energy += frequency * frequency;
      entropy -= ( frequency > 0.0001 ) ? frequency *std::log(frequency) / log2:0;
      correlation += ( ( a - pixelMean ) * ( b - pixelMean ) * frequency ) / pixelVarianceSquared;
      inverseDifferenceMoment += frequency / ( 1.0 + ( a - b ) * ( a - b ) );
      inertia += ( a - b ) * ( a - b ) * frequency;
      clusterShade += std::pow( ( a - pixelMean ) + ( b - pixelMean ), 3 )  * frequency;
      clusterProminence += std::pow( ( a - pixelMean ) + ( b - pixelMean ), 4 ) * frequency;
      haralickCorrelation += a * b * frequency;
        }
      }

    haralickCorrelation = ( haralickCorrelation - marginalMean * marginalMean ) / marginalDevSquared;

    outputPixel[0] = energy;
    outputPixel[1] = entropy;
    outputPixel[2] = correlation;
    outputPixel[3] = inverseDifferenceMoment;
    outputPixel[4] = inertia;
    outputPixel[5] = clusterShade;
    outputPixel[6] = clusterProminence;
    outputPixel[7] = haralickCorrelation;
}

template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::ComputeMeansAndVariances(const vnl_matrix<unsigned int> &hist,
                           const unsigned int totalNumberOfFreq,
                           double & pixelMean,
                           double & marginalMean,
                           double & marginalDevSquared,
                           double & pixelVariance)
{
  // This function takes two passes through the histogram and two passes through
  // an array of the same length as a histogram axis. This could probably be
  // cleverly compressed to one pass, but it's not clear that that's necessary.

  // Initialize everything
  double *marginalSums = new double[m_NumberOfBinsPerAxis];

  for ( double *ms_It = marginalSums;
        ms_It < marginalSums + m_NumberOfBinsPerAxis; ms_It++ )
    {
    *ms_It = 0;
    }
  pixelMean = 0;

  // Ok, now do the first pass through the histogram to get the marginal sums
  // and compute the pixel mean
  for(unsigned int a = 0; a < m_NumberOfBinsPerAxis; a++)
    {
    for(unsigned int b = 0; b < m_NumberOfBinsPerAxis; b++)
      {
      float frequency = hist[a][b] / (float)totalNumberOfFreq;
      pixelMean += a * frequency;
      marginalSums[a] += frequency;
      }
    }

  /*  Now get the mean and deviaton of the marginal sums.
      Compute incremental mean and SD, a la Knuth, "The  Art of Computer
      Programming, Volume 2: Seminumerical Algorithms",  section 4.2.2.
      Compute mean and standard deviation using the recurrence relation:
      M(1) = x(1), M(k) = M(k-1) + (x(k) - M(k-1) ) / k
      S(1) = 0, S(k) = S(k-1) + (x(k) - M(k-1)) * (x(k) - M(k))
      for 2 <= k <= n, then
      sigma = std::sqrt(S(n) / n) (or divide by n-1 for sample SD instead of
      population SD).
  */
  marginalMean = marginalSums[0];
  marginalDevSquared = 0;
  for ( unsigned int arrayIndex = 1; arrayIndex < m_NumberOfBinsPerAxis; arrayIndex++ )
    {
    int    k = arrayIndex + 1;
    double M_k_minus_1 = marginalMean;
    double S_k_minus_1 = marginalDevSquared;
    double x_k = marginalSums[arrayIndex];

    double M_k = M_k_minus_1 + ( x_k - M_k_minus_1 ) / k;
    double S_k = S_k_minus_1 + ( x_k - M_k_minus_1 ) * ( x_k - M_k );

    marginalMean = M_k;
    marginalDevSquared = S_k;
    }
  marginalDevSquared = marginalDevSquared / m_NumberOfBinsPerAxis;

  // OK, now compute the pixel variances.
  pixelVariance = 0;
  for(unsigned int a = 0; a < m_NumberOfBinsPerAxis; a++)
    {
    for(unsigned int b = 0; b < m_NumberOfBinsPerAxis; b++)
      {
      float frequency = hist[a][b] / (float)totalNumberOfFreq;
      pixelVariance += ( a - pixelMean ) * ( a - pixelMean ) * (frequency);
      }
    }

  delete[] marginalSums;
}

template<typename TInputImage, typename TOutputImage>
void
CoocurrenceTextureFeaturesImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{

  Superclass::PrintSelf( os, indent );

  itkPrintSelfObjectMacro( DigitalizedInputImage );

  os << indent << "NeighborhoodRadius"
    << static_cast< typename NumericTraits<
    NeighborhoodRadiusType >::PrintType >( m_NeighborhoodRadius ) << std::endl;

  itkPrintSelfObjectMacro( Offsets );

  os << indent << "NumberOfBinsPerAxis" << m_NumberOfBinsPerAxis << std::endl;
  os << indent << "Min"
    << static_cast< typename NumericTraits< PixelType >::PrintType >( m_Min )
    << std::endl;
  os << indent << "Max"
    << static_cast< typename NumericTraits< PixelType >::PrintType >( m_Max )
    << std::endl;
  os << indent << "InsidePixelValue"
    << static_cast< typename NumericTraits< PixelType >::PrintType >(
    m_InsidePixelValue ) << std::endl;
  os << indent << "Spacing"
    << static_cast< typename NumericTraits<
    typename TInputImage::SpacingType >::PrintType >( m_Spacing ) << std::endl;
  os << indent << "Normalize" << m_Normalize << std::endl;
}
} // end of namespace Statistics
} // end of namespace itk

#endif
