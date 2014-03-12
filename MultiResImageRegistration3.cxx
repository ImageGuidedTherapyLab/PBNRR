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

//
// This example is the 3D version of the 2D example in MultiResImageRegistration1.cxx
//

#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkTranslationTransform.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCheckerBoardImageFilter.h"


#include "itkCommand.h"


template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );

protected:
  RegistrationInterfaceCommand() {};

public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
  typedef   OptimizerType *                            OptimizerPointer;

  void Execute(itk::Object * object, const itk::EventObject & event)
    {
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
      return;
      }
    RegistrationPointer registration =
      dynamic_cast<RegistrationPointer>( object );

    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >(
      registration->GetModifiableOptimizer() );

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : "
              << registration->GetCurrentLevel()  << std::endl;
    std::cout << std::endl;

    if ( registration->GetCurrentLevel() == 0 )
      {
      optimizer->SetMaximumStepLength( 16.00 );
      optimizer->SetMinimumStepLength( 0.01 );
      }
    else
      {
      optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / 4.0 );
      optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 10.0 );
      }
    }

  void Execute(const itk::Object * , const itk::EventObject & )
    { return; }
};


//  The following section of code implements an observer
//  that will monitor the evolution of the registration process.
//
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:
  typedef   itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef   const OptimizerType *                    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer =
        dynamic_cast< OptimizerPointer >( object );
      if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};


int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << " outputImagefile ";
    std::cerr << " [maxNumberOfIterations] [resolutionLevels] " << std::endl;
    std::cerr << " [translationScale] [relaxationFactor] " << std::endl;
    std::cerr << " [numberOfBins] [numberOfSamples ] " << std::endl;
    std::cerr << " [useExplicitPDFderivatives ] " << std::endl;
    std::cerr << " [checkerBoardBefore] [checkerBoardAfter]";
    return EXIT_FAILURE;
    }
  enum { FIXED_IMG = 1, MOVING_IMG, OUTPUT_IMG, MAX_ITER, NUM_LEVELS, TRANS_SCALE, RELAX_FACT, NUM_BINS, NUM_SAMPLE, EXP_DERIV, CHECKERBEFORE, CHECKERAFTER};

  const    unsigned int    Dimension = 3;
  typedef  unsigned short  PixelType;

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  typedef   float                                    InternalPixelType;
  typedef itk::Image< InternalPixelType, Dimension > InternalImageType;

  //typedef itk::TranslationTransform< double, Dimension > TransformType;
  typedef itk::AffineTransform< double, Dimension > TransformType;
  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
  typedef itk::LinearInterpolateImageFunction<
                                    InternalImageType,
                                    double             > InterpolatorType;
  typedef itk::MattesMutualInformationImageToImageMetric<
                                    InternalImageType,
                                    InternalImageType >   MetricType;
  typedef itk::MultiResolutionImageRegistrationMethod<
                                    InternalImageType,
                                    InternalImageType >   RegistrationType;

  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   FixedImagePyramidType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   MovingImagePyramidType;


  //  All the components are instantiated using their \code{New()} method
  //  and connected to the registration object as in previous example.
  //
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  MetricType::Pointer         metric        = MetricType::New();

  FixedImagePyramidType::Pointer fixedImagePyramid =
      FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingImagePyramid =
      MovingImagePyramidType::New();

  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetInterpolator(  interpolator  );
  registration->SetMetric( metric  );
  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );


  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[FIXED_IMG ] );
  movingImageReader->SetFileName( argv[MOVING_IMG] );


  typedef itk::CastImageFilter<
                        FixedImageType, InternalImageType > FixedCastFilterType;
  typedef itk::CastImageFilter<
                        MovingImageType, InternalImageType > MovingCastFilterType;

  FixedCastFilterType::Pointer fixedCaster   = FixedCastFilterType::New();
  MovingCastFilterType::Pointer movingCaster = MovingCastFilterType::New();

  fixedCaster->SetInput(  fixedImageReader->GetOutput() );
  movingCaster->SetInput( movingImageReader->GetOutput() );

  registration->SetFixedImage(    fixedCaster->GetOutput()    );
  registration->SetMovingImage(   movingCaster->GetOutput()   );


  fixedCaster->Update();

  registration->SetFixedImageRegion(
       fixedCaster->GetOutput()->GetBufferedRegion() );


  typedef RegistrationType::ParametersType ParametersType;
  //ParametersType initialParameters( transform->GetNumberOfParameters() );

  //initialParameters[0] = 0.0;  // Initial offset in mm along X
  //initialParameters[1] = 0.0;  // Initial offset in mm along Y
  //initialParameters[2] = 0.0;  // Initial offset in mm along Z

  //registration->SetInitialTransformParameters( initialParameters );
  registration->SetInitialTransformParameters(
                                 transform->GetParameters() );

  //  Software Guide : BeginLatex
  //
  //  Keeping in mind that the scale of units in scaling, rotation and
  //  translation are quite different, we take advantage of the scaling
  //  functionality provided by the optimizers. We know that the first $N
  //  \times N$ elements of the parameters array correspond to the rotation
  //  matrix factor, and the last $N$ are the components of the translation to
  //  be applied after multiplication with the matrix is performed.
  //
  //  Software Guide : EndLatex


  double translationScale = 1.0 / 1000.0;
  if( argc > TRANS_SCALE )
    {
    translationScale = atof( argv[TRANS_SCALE] );
    }


  // Software Guide : BeginCodeSnippet
  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

  optimizerScales[0] =  1.0;
  optimizerScales[1] =  1.0;
  optimizerScales[2] =  1.0;
  optimizerScales[3] =  1.0;
  optimizerScales[4] =  1.0;
  optimizerScales[5] =  1.0;
  optimizerScales[6] =  1.0;
  optimizerScales[7] =  1.0;
  optimizerScales[8] =  1.0;
  optimizerScales[9]  =  translationScale;
  optimizerScales[10] =  translationScale;
  optimizerScales[11] =  translationScale;

  optimizer->SetScales( optimizerScales );
  // Software Guide : EndCodeSnippet

  metric->SetNumberOfHistogramBins( 128 );
  metric->SetNumberOfSpatialSamples( 50000 );

  if( argc > NUM_BINS )
    {
    // optionally, override the values with numbers taken from the command line arguments.
    metric->SetNumberOfHistogramBins( atoi( argv[NUM_BINS] ) );
    }

  if( argc > NUM_SAMPLE )
    {
    // optionally, override the values with numbers taken from the command line arguments.
    metric->SetNumberOfSpatialSamples( atoi( argv[NUM_SAMPLE] ) );
    }

  metric->ReinitializeSeed( 76926294 );


  if( argc > EXP_DERIV )
    {
    // Define whether to calculate the metric derivative by explicitly
    // computing the derivatives of the joint PDF with respect to the Transform
    // parameters, or doing it by progressively accumulating contributions from
    // each bin in the joint PDF.
    metric->SetUseExplicitPDFDerivatives( atoi( argv[EXP_DERIV] ) );
    }


  optimizer->SetNumberOfIterations( 200 );
  optimizer->SetRelaxationFactor( 0.9 );

  if( argc > MAX_ITER )
    {
    // optionally, override the values with numbers taken from the command line arguments.
    optimizer->SetNumberOfIterations( atoi( argv[MAX_ITER] ) );
    }

  if( argc > RELAX_FACT )
    {
    // optionally, override the values with numbers taken from the command line arguments.
    optimizer->SetRelaxationFactor( atof( argv[RELAX_FACT] ) );
    }

  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
  CommandType::Pointer command = CommandType::New();
  registration->AddObserver( itk::IterationEvent(), command );


  registration->SetNumberOfLevels( 3 );

  try
    {
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
    }

  ParametersType finalParameters = registration->GetLastTransformParameters();

  //double TranslationAlongX = finalParameters[0];
  //double TranslationAlongY = finalParameters[1];
  //double TranslationAlongZ = finalParameters[2];

  unsigned int numberOfIterations = optimizer->GetCurrentIteration();

  double bestValue = optimizer->GetValue();


  // Print out results
  //
  std::cout << "Result = " << finalParameters << std::endl;
  //std::cout << " Translation X = " << TranslationAlongX  << std::endl;
  //std::cout << " Translation Y = " << TranslationAlongY  << std::endl;
  //std::cout << " Translation Z = " << TranslationAlongZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetParameters( finalParameters );
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( finalTransform );
  resample->SetInput( movingImageReader->GetOutput() );

  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

  PixelType backgroundGrayLevel = 0;

  resample->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( backgroundGrayLevel );


  typedef  short          OutputPixelType;

  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter<
                        FixedImageType,
                        OutputImageType > CastFilterType;

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;


  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();


  writer->SetFileName( argv[OUTPUT_IMG] );


  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  //
  // Generate checkerboards before and after registration
  //
  typedef itk::CheckerBoardImageFilter< FixedImageType > CheckerBoardFilterType;

  CheckerBoardFilterType::Pointer checker = CheckerBoardFilterType::New();

  checker->SetInput1( fixedImage );
  checker->SetInput2( resample->GetOutput() );

  caster->SetInput( checker->GetOutput() );
  writer->SetInput( caster->GetOutput()   );

  resample->SetDefaultPixelValue( 0 );

  // Before registration
  TransformType::Pointer identityTransform = TransformType::New();
  identityTransform->SetIdentity();
  resample->SetTransform( identityTransform );

  if( argc > CHECKERBEFORE )
    {
    writer->SetFileName( argv[CHECKERBEFORE] );
    writer->Update();
    }


  // After registration
  resample->SetTransform( finalTransform );
  if( argc > CHECKERAFTER )
    {
    writer->SetFileName( argv[CHECKERAFTER] );
    writer->Update();
    }

  return EXIT_SUCCESS;
}
