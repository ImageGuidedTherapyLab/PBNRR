//=========================================================================
// *
// *  Copyright Insight Software Consortium
// *
// *  Licensed under the Apache License, Version 2.0 (the "License");
// *  you may not use this file except in compliance with the License.
// *  You may obtain a copy of the License at
// *
// *         http://www.apache.org/licenses/LICENSE-2.0.txt
// *
// *  Unless required by applicable law or agreed to in writing, software
// *  distributed under the License is distributed on an "AS IS" BASIS,
// *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// *  See the License for the specific language governing permissions and
// *  limitations under the License.
// *
// *=========================================================================*/

#include <iostream>

#include "itkVTKTetrahedralMeshReader.h"
#include "itkImage.h"
//#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkPhysicsBasedNonRigidRegistrationMethod.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <sys/time.h>

const unsigned int ImageDimension = 3;
typedef float                                                     InputPixelType;
typedef short                                                   MaskPixelType;

typedef itk::Matrix< double , ImageDimension, ImageDimension >        PointSetPixelType;
typedef itk::PointSet< PointSetPixelType, ImageDimension >          PointSetType;

typedef itk::Image< InputPixelType,  ImageDimension >             InputImageType;
typedef itk::Image< MaskPixelType,  ImageDimension >              MaskImageType;
typedef itk::ImageFileReader< InputImageType >                    ImageReaderType;
typedef itk::Image< itk::Vector< float, ImageDimension >, ImageDimension >  DeformationFieldType;
typedef itk::Mesh< float, ImageDimension >                      MeshType;

typedef itk::fem::PhysicsBasedNonRigidRegistrationMethod<InputImageType, InputImageType, InputImageType, MeshType, DeformationFieldType>  PBNRRFilterType;
typedef itk::MaskFeaturePointSelectionFilter< InputImageType,InputImageType, PointSetType >  FeatureSelectionFilterType;
typedef itk::BlockMatchingImageFilter< InputImageType >   BlockMatchingFilterType;



int main(int argc, char *argv[] )
{

  if ( argc < 6)
  {
    std::cerr << "Five arguments are required :"<< std::endl;
    std::cerr <<" FixedImage, MovingImage, MaskImage, Mesh, WarpedImage" << std::endl;
    return EXIT_FAILURE;
  }

  enum { FIXED_IMG = 1, MOVING_IMG, MASK_IMG, MESH, WARPED_IMG };

  // time
  timeval tim;
  gettimeofday(&tim, NULL);
  double startTotalTime = tim.tv_sec + (tim.tv_usec/1000000.0);

  // read fixed image
  ImageReaderType::Pointer readerFixed = ImageReaderType::New();
  readerFixed->SetFileName( argv[FIXED_IMG] );

  // read moving image
  ImageReaderType::Pointer readerMoving = ImageReaderType::New();
  readerMoving->SetFileName( argv[MOVING_IMG] );

  // read mask image
  ImageReaderType::Pointer readerMask = ImageReaderType::New();
  readerMask->SetFileName( argv[MASK_IMG] );

  // read mesh
  typedef itk::VTKTetrahedralMeshReader< MeshType >   MeshReaderType;
  MeshReaderType::Pointer readerMesh = MeshReaderType::New();
  readerMesh->SetFileName( argv[MESH] );

  // Update the readers
  try
  {
    readerFixed->Update();
    readerMoving->Update();
    readerMask->Update();
    readerMesh->Update();
  }
  catch( itk::ExceptionObject & e )
  {
    std::cerr << "Error while reading inputs: " << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
  }

  // input prints
  std::cout << "*********************************************" << std::endl;
  std::cout << "Moving Image : " << argv[MOVING_IMG] << std::endl;
  std::cout << "Fixed Image : " << argv[FIXED_IMG] << std::endl;
  std::cout << "Mask Image : " << argv[MASK_IMG] << std::endl;
  std::cout << "Mesh : " << argv[MESH] << std::endl;
  std::cout << "Warped Image : " << argv[WARPED_IMG] << std::endl;
  std::cout << "*********************************************" << std::endl;
  std::cout << " " << std::endl;

  // Create the NRR filter and set the input  
  PBNRRFilterType::Pointer filter = PBNRRFilterType::New();
  filter->SetFixedImage( readerFixed->GetOutput() );
  filter->SetMovingImage( readerMoving->GetOutput() );
  filter->SetMaskImage( readerMask->GetOutput() );
  filter->SetMesh( readerMesh->GetOutput() );

  itk::Size< ImageDimension > BlockRadious;
  BlockRadious.Fill(1);
  filter->SetBlockRadius(BlockRadious);

  itk::Size< ImageDimension > SearchRadious;
  SearchRadious.Fill(5);
  filter->SetSearchRadius(SearchRadious);

  filter->SetApproximationSteps(10);
  filter->SetOutlierRejectionSteps(10);
  filter->SetSelectFraction( 0.05 );

  std::cout << "Filter: " << filter << std::endl;

  // Update the NRR filter
  try
  {
    filter->Update();
  }
  catch( itk::ExceptionObject & e )
  {
    std::cerr << "Error during filter->Update(): " << e << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Done." << std::endl;

  // time
  gettimeofday(&tim, NULL);
  double startTime = tim.tv_sec + (tim.tv_usec/1000000.0);

  // // Create - Write ITK deformed image
  // InputImageType::Pointer deformedImage;
  // filter->CreateDeformedImage(deformedImage);
  // 
  // std::cout << "Save Deformed Image at  : " << argv[WARPED_IMG] << std::endl;
  // typedef itk::ImageFileWriter<InputImageType> WriterType;
  // typename WriterType::Pointer deformedImageWriter = WriterType::New();
  // deformedImageWriter->SetFileName(argv[WARPED_IMG]);
  // deformedImageWriter->SetInput(deformedImage);
  // deformedImageWriter->Update();

  gettimeofday(&tim, NULL);
  double endTime  = tim.tv_sec + (tim.tv_usec/1000000.0);
  std::cout << "" << std::endl;
  std::cout << "Time for creating the Deformed Image : " << endTime - startTime <<  " sec" << std::endl ;
  std::cout << "" << std::endl;

  // Print time
  gettimeofday(&tim, NULL);
  double endTotalTime  = tim.tv_sec + (tim.tv_usec/1000000.0);
  std::cout << "" << std::endl ;
  std::cout << "Total running time : " << endTotalTime - startTotalTime <<  " sec" << std::endl ;
  std::cout << "" << std::endl ;

  return EXIT_SUCCESS;

}
