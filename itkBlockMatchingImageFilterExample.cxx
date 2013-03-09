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

#include <iostream>
#include "itkIndex.h"
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLineIterator.h"
#include "itkMultiThreader.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkMaskFeaturePointSelectionFilter.h"
#include "itkBlockMatchingVaryImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkResampleImageFilter.h"
#include "GetPot.h"


int main( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl<< std::endl;
    std::cerr << " itkitkBlockMatchingImageFilterExample IniFile" << std::endl << std::endl;
    std::cerr << "   step one: find all feature points on a moving image" << std::endl;
    std::cerr << "             within a variance threshold              " << std::endl;
    std::cerr << "   step two: within search window on fixed image      " << std::endl;
    std::cerr << "             find block radius window center location " << std::endl;
    std::cerr << "             with the highest similarity to the moving" << std::endl;
    std::cerr << "             image feature point" << std::endl << std::endl;
    std::cerr << "   assumption: the search window on the fixed image   " << std::endl;
    std::cerr << "               is the feature point location in       " << std::endl;
    std::cerr << "               physical space +- the search radius    " << std::endl;
    std::cerr << "               ie, all point match pairs should be    " << std::endl;
    std::cerr << "               within the search radius neighborhood  " << std::endl;
    std::cerr << "               defined the the physical space location" << std::endl;
    std::cerr << "               of the feature point +- the search radius" << std::endl;
    std::cerr << "   assumption: the fixed and moving image should have " << std::endl;
    std::cerr << "               the same dimensions    " << std::endl<< std::endl;
    std::cerr << "   ...ie, rigid registration and resampling are probably needed" << std::endl << std::endl;
    return EXIT_FAILURE;
    }

  GetPot controlfile( argv[1] );

  const double selectFraction = controlfile("featureselection/selectionfraction",0.01);

  typedef unsigned char                  InputPixelType;
  typedef double                         NumpyArrayType;
  static const unsigned int Dimension = 3;
  typedef itk::Image< InputPixelType,  Dimension >  InputImageType;

  //typedef itk::RGBPixel<InputPixelType>  OutputPixelType;
  //typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;

  // set anisotropic match and search windows
  // ini file should of the form
  // [blockmatch]
  // blockradius  = '2 2 3'
  // searchradius = '5 5 7'
  // 
  // Parameters used for FS and BM
  typedef InputImageType::SizeType RadiusType;
  RadiusType blockRadius, searchRadius;
  for (int icoord = 0; icoord < Dimension; icoord++)
   {
    blockRadius.SetElement(  icoord, 
                             controlfile("blockmatch/blockradius" , 3, icoord) );
    searchRadius.SetElement( icoord,
                             controlfile("blockmatch/searchradius", 7, icoord) );
   }


  typedef itk::ImageFileReader< InputImageType >  ReaderType;

  // read in moving image first and identify the feature points
  ReaderType::Pointer readerMoving = ReaderType::New();
  readerMoving->SetFileName( controlfile("image/moving","./MovingNotFound") );
  try
    {
    readerMoving->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Error in reading the input image: " << e << std::endl;
    return EXIT_FAILURE;
    }

  // Reduce region of interest by SEARCH_RADIUS
  typedef itk::RegionOfInterestImageFilter< InputImageType, InputImageType >  RegionOfInterestFilterType;

  RegionOfInterestFilterType::Pointer regionOfInterestFilter = RegionOfInterestFilterType::New();

  regionOfInterestFilter->SetInput( readerMoving->GetOutput() );

  RegionOfInterestFilterType::RegionType regionOfInterest = readerMoving->GetOutput()->GetLargestPossibleRegion();

  // set the starting point of the ROI
  RegionOfInterestFilterType::RegionType::IndexType regionOfInterestIndex = regionOfInterest.GetIndex();
  regionOfInterestIndex += searchRadius;
  regionOfInterest.SetIndex( regionOfInterestIndex );

  // the ROI size is the full size minus the radius on the beging
  //   and minus the radius on the end 
  RegionOfInterestFilterType::RegionType::SizeType regionOfInterestSize = regionOfInterest.GetSize();
  regionOfInterestSize -= searchRadius + searchRadius;
  regionOfInterest.SetSize( regionOfInterestSize );

  regionOfInterestFilter->SetRegionOfInterest( regionOfInterest );
  regionOfInterestFilter->Update();

  regionOfInterestFilter->DebugOn();
  if ( regionOfInterestFilter->GetDebug() )
   {
     std::ostringstream ROIFileName;
     ROIFileName  << controlfile("image/output","./Output") << "ROI.mha" ;
     //Set up the writer
     typedef itk::ImageFileWriter< InputImageType >  WriterType;
     WriterType::Pointer writer = WriterType::New();

     writer->SetFileName( ROIFileName.str() );
     writer->SetInput( regionOfInterestFilter->GetOutput() );
     try
       {
       writer->Update();
       }
     catch( itk::ExceptionObject & e )
       {
       std::cerr << "Error in writing the output image:" << e << std::endl;
       return EXIT_FAILURE;
       }
   }

  // typedefs
  typedef itk::MaskFeaturePointSelectionFilter< InputImageType >  FeatureSelectionFilterType;
  typedef FeatureSelectionFilterType::FeaturePointsType           PointSetType;

  typedef FeatureSelectionFilterType::PointType       PointType;
  typedef FeatureSelectionFilterType::InputImageType  ImageType;

  // Feature Selection
  FeatureSelectionFilterType::Pointer featureSelectionFilter = FeatureSelectionFilterType::New();

  featureSelectionFilter->SetInput( regionOfInterestFilter->GetOutput() );
  featureSelectionFilter->SetSelectFraction( selectFraction );
  featureSelectionFilter->SetBlockRadius( blockRadius );
  featureSelectionFilter->ComputeStructureTensorsOff();
  featureSelectionFilter->SetNonConnectivity((unsigned int)controlfile("featureselection/nonconnectivity",0));

  // check if user input mask image
  std::string MovingMaskFileName( controlfile("image/movingmask","MovingMaskNotFound") ) ;
  if (MovingMaskFileName.find("MovingMaskNotFound")!=std::string::npos)
   {
     std::cout << "No Moving Image Mask input... " << std::endl;
     std::cout << "No Moving Image Mask input... " << std::endl;
   }
  else
   {
    ReaderType::Pointer readerMovingMask = ReaderType::New();
    readerMovingMask->SetFileName( MovingMaskFileName.c_str() );
    try
      {
        readerMovingMask->Update();
      }
    catch( itk::ExceptionObject & e )
      {
        std::cerr << "Error in reading the moving image mask image: " <<  MovingMaskFileName << e << std::endl;
      }
    // set mask image
    featureSelectionFilter->SetMaskImage( readerMovingMask->GetOutput() );
   }
  // update and write selection points
  featureSelectionFilter->Update();
  // open file to write feature points
  std::ofstream      FeatureFile; 
  std::ostringstream FeatureFileName;
  FeatureFileName       << controlfile("image/output","./Output") <<  "Feature.txt"      ;
  FeatureFile.open      ( FeatureFileName.str().c_str()       );

  typedef itk::Point<float, Dimension >      ITKFloatPointType;
  typedef std::vector< ITKFloatPointType > STLFeaturePointType;
  const STLFeaturePointType &stlfeaturepoints = 
        featureSelectionFilter->GetOutput()->GetPoints()->CastToSTLConstContainer();
  const int NumFeaturePoints = stlfeaturepoints.size(); 
  // pass back to python as 1d array, reshape in numpy...
  std::vector< NumpyArrayType > StlFeaturePointsPassToCython;

  // write points
  for (STLFeaturePointType::const_iterator ptIt = stlfeaturepoints.begin(); 
                                           ptIt!= stlfeaturepoints.end()  ; ++ptIt) 
    {
    const ITKFloatPointType &featurepoint = *ptIt;
    for (int icoord = 0; icoord < Dimension; icoord++)
      StlFeaturePointsPassToCython.push_back(  featurepoint[icoord] );
    FeatureFile << featurepoint[0]<<" "
                << featurepoint[1]<<" "
                << featurepoint[2]<< std::endl;
    }
  // close file
  FeatureFile.close();
  std::cout << "Wrote "<< NumFeaturePoints << " Feature Points..." << std::endl;

  //Set up the reader
  ReaderType::Pointer readerFixed = ReaderType::New();
  readerFixed->SetFileName( controlfile("image/fixed","./FixedNotFound") );
  try
    {
    readerFixed->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Error in reading the input image: " << e << std::endl;
    return EXIT_FAILURE;
    }


  typedef itk::BlockMatchingVaryImageFilter< InputImageType >  BlockMatchingFilterType;
  BlockMatchingFilterType::Pointer blockMatchingFilter = BlockMatchingFilterType::New();

  // inputs (all required)
  blockMatchingFilter->SetFixedImage(  readerFixed->GetOutput() );
  blockMatchingFilter->SetMovingImage( readerMoving->GetOutput() );
  blockMatchingFilter->SetFeaturePoints( featureSelectionFilter->GetOutput() );

  // parameters (all optional)
  blockMatchingFilter->SetNumberOfThreads( controlfile("exec/threads" , 
       (int)  blockMatchingFilter->GetMultiThreader()->GetGlobalDefaultNumberOfThreads()
                                         )  );
  //blockMatchingFilter->DebugOn( );
  blockMatchingFilter->SetBlockRadius( blockRadius );
  blockMatchingFilter->SetSearchRadius( searchRadius );

  std::cout << "Block matching: " << blockMatchingFilter << std::endl;
  try
    {
    blockMatchingFilter->Update();
    }
  catch ( itk::ExceptionObject &err )
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Exercise the following methods
  BlockMatchingFilterType::DisplacementsType * displacements = blockMatchingFilter->GetDisplacements();
  if( displacements == NULL )
    {
    std::cerr << "GetDisplacements() failed." << std::endl;
    return EXIT_FAILURE;
    }
  BlockMatchingFilterType::SimilaritiesType * similarities = blockMatchingFilter->GetSimilarities();
  if( similarities == NULL )
    {
    std::cerr << "GetSimilarities() failed." << std::endl;
    return EXIT_FAILURE;
    }

  // // create RGB copy of input image
  // typedef itk::ScalarToRGBColormapImageFilter< InputImageType, OutputImageType >  RGBFilterType;
  // RGBFilterType::Pointer colormapImageFilter = RGBFilterType::New();

  // colormapImageFilter->SetColormap( RGBFilterType::Grey );
  // colormapImageFilter->SetInput( readerMoving->GetOutput() );
  // try
  //   {
  //   colormapImageFilter->Update();
  //   }
  // catch ( itk::ExceptionObject &err )
  //   {
  //   std::cerr << err << std::endl;
  //   return EXIT_FAILURE;
  //   }

  // OutputImageType::Pointer outputImage = colormapImageFilter->GetOutput();

  // Highlight the feature points identified in the output image

  // // define colors
  // OutputPixelType red;
  // red.SetRed( 255 );
  // red.SetGreen( 0 );
  // red.SetBlue( 0 );

  // OutputPixelType green;
  // green.SetRed( 0 );
  // green.SetGreen( 255 );
  // green.SetBlue( 0 );

  // OutputPixelType blue;
  // blue.SetRed( 0 );
  // blue.SetGreen( 0 );
  // blue.SetBlue( 255 );

  // OutputImageType::IndexType index;
  // while ( pointItr != pointEnd )
  //   {
  //   if ( outputImage->TransformPhysicalPointToIndex(pointItr.Value(), index) )
  //     {
  //     OutputImageType::IndexType displ;
  //     outputImage->TransformPhysicalPointToIndex( pointItr.Value() + displItr.Value(), displ );

  //     // draw line between old and new location of a point in blue
  //     itk::LineIterator< OutputImageType > lineIter( outputImage, index, displ );
  //     for ( lineIter.GoToBegin(); !lineIter.IsAtEnd(); ++lineIter )
  //       {
  //       lineIter.Set( blue );
  //       }

  //     // mark old location of a point in green
  //     outputImage->SetPixel(index, green);

  //     // mark new location of a point in red
  //     outputImage->SetPixel(displ, red);
  //     }
  //   pointItr++;
  //   displItr++;
  //   }

  // write txt files of block matching data
  //   source, target, displacements, similarity
  std::ofstream      DisplacementFile,    SimilarityFile;
  std::ostringstream DisplacementFileName,SimilarityFileName;

  DisplacementFileName << controlfile("image/output","./Output") <<  "Displacement.txt";
  SimilarityFileName   << controlfile("image/output","./Output") <<  "Similarity.txt"  ;

  // similarities iterator
  const std::vector< NumpyArrayType > &stlsimilarities = similarities->GetPointData()->CastToSTLConstContainer();
  const int NumSimilarityPoints = stlsimilarities.size(); 

  // displacement iterator
  typedef itk::Vector<float, Dimension >            ITKFloatVectorType;
  typedef std::vector< ITKFloatVectorType > STLDisplacementVectorType;
  const STLDisplacementVectorType &stldisplacements =  displacements->GetPointData()->CastToSTLConstContainer();
  const int NumDisplacementVec = stldisplacements.size(); 

  // error check
  if( NumDisplacementVec  != NumFeaturePoints 
             or
      NumFeaturePoints    != NumSimilarityPoints 
             or
      NumSimilarityPoints != NumDisplacementVec  
    )
    {
    std::cerr << "Number of points do not match !!! "           << std::endl;
    std::cerr << "NumFeaturePoints:    " << NumFeaturePoints    << std::endl;
    std::cerr << "NumSimilarityPoints: " << NumSimilarityPoints << std::endl;
    std::cerr << "NumDisplacementVec:  " << NumDisplacementVec  << std::endl;
    return EXIT_FAILURE;
    }

  // open file and write similarity
  SimilarityFile.open  ( SimilarityFileName.str().c_str()   );
  for (std::vector< NumpyArrayType >::const_iterator similItr = stlsimilarities.begin(); 
                                             similItr!= stlsimilarities.end()  ; ++similItr) 
    {
    SimilarityFile   << *similItr << std::endl;
    }
  SimilarityFile.close(); 


  // pass back to python as 1d array, reshape in numpy...
  std::vector< NumpyArrayType > StlDisplacementsPassToCython;
  // open file and write displacement
  DisplacementFile.open( DisplacementFileName.str().c_str() );
  for (STLDisplacementVectorType::const_iterator dplIt = stldisplacements.begin(); 
                                                 dplIt!= stldisplacements.end()  ; ++dplIt) 
    {
    const ITKFloatVectorType &displacementvalue = *dplIt;
    for (int icoord = 0; icoord < Dimension; icoord++)
      StlDisplacementsPassToCython.push_back(  displacementvalue[icoord] );
    DisplacementFile << displacementvalue[0]<<" "
                     << displacementvalue[1]<<" "
                     << displacementvalue[2]<< std::endl;
    }
  DisplacementFile.close();

  // //Set up the writer
  // typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  // WriterType::Pointer writer = WriterType::New();

  // writer->SetFileName( argv[3] );
  // writer->SetInput( outputImage );
  // try
  //   {
  //   writer->Update();
  //   }
  // catch( itk::ExceptionObject & e )
  //   {
  //   std::cerr << "Error in writing the output image:" << e << std::endl;
  //   return EXIT_FAILURE;
  //   }

  return EXIT_SUCCESS;
}
