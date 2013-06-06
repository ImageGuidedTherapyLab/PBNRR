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

#include <algorithm>

#include "itkCommand.h"
#include "itkImageRegistrationMethod.h"
#include "itkExhaustiveOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkTranslationTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "vnl/vnl_math.h"

/**
 *  The objectif function is the quadratic form:
 *
 *  1/2 x^T A x - b^T x
 *
 *  Where A is a matrix and b is a vector
 *  The system in this example is:
 *
 *     | 3  2 ||x|   | 2|   |0|
 *     | 2  6 ||y| + |-8| = |0|
 *
 *
 *   the solution is the vector | 2 -2 |
 *
 * \class RSGCostFunction
 */
class RSGCostFunction : public itk::SingleValuedCostFunction
{
public:

  typedef RSGCostFunction                 Self;
  typedef itk::SingleValuedCostFunction   Superclass;
  typedef itk::SmartPointer<Self>         Pointer;
  typedef itk::SmartPointer<const Self>   ConstPointer;
  itkNewMacro( Self );

  enum { SpaceDimension=2 };

  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef Superclass::MeasureType         MeasureType;

  RSGCostFunction()
  {
  }


  MeasureType  GetValue( const ParametersType & parameters ) const
  {

    double x = parameters[0];
    double y = parameters[1];

    std::cout << "GetValue( ";
    std::cout << x << " ";
    std::cout << y << ") = ";

    MeasureType measure = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y;

    std::cout << measure << std::endl;

    return measure;

  }

  void GetDerivative( const ParametersType & parameters,
                            DerivativeType  & derivative ) const
  {

    double x = parameters[0];
    double y = parameters[1];

    std::cout << "GetDerivative( ";
    std::cout << x << " ";
    std::cout << y << ") = ";

    derivative = DerivativeType( SpaceDimension );
    derivative[0] = 3 * x + 2 * y -2;
    derivative[1] = 2 * x + 6 * y +8;

  }


  unsigned int GetNumberOfParameters(void) const
    {
    return SpaceDimension;
    }
};

class IndexObserver : public itk::Command
{
public:
  typedef IndexObserver              Self;
  typedef itk::Command               Superclass;
  typedef itk::SmartPointer < Self > Pointer;

  itkNewMacro ( IndexObserver );

  virtual void  Execute ( const itk::Object *caller, const itk::EventObject &)
  {
    typedef itk::ExhaustiveOptimizer OptimizerType;
    const OptimizerType *optimizer = dynamic_cast < const OptimizerType * > ( caller );

    if ( 0 != optimizer )
    {
      OptimizerType::ParametersType currentIndex = optimizer->GetCurrentIndex ();

      if ( currentIndex.GetSize () == 2 )
      {
        std::cout << " @ index = " << currentIndex << std::endl;
        // Casting is safe here since the indices are always integer values (but there are stored in doubles):
        unsigned long idx = static_cast < unsigned long > ( currentIndex [ 0 ] + 21 * currentIndex [ 1 ] );
        m_VisitedIndices.push_back ( idx );
      }
    }
  }

  virtual void  Execute (itk::Object *caller, const itk::EventObject &event)
  {
    Execute ( static_cast < const itk::Object * > ( caller ), event );
  }

  std::vector < unsigned long > m_VisitedIndices;
};

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile output" << std::endl;
    return EXIT_FAILURE;
    }
  enum { FIXED_IMG = 1, MOVING_IMG, OUTPUT_IMG};

  const    unsigned int    Dimension = 3;
  typedef   float                                    InternalPixelType;
  typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
  typedef itk::ImageFileReader< InternalImageType > FixedImageReaderType;
  typedef itk::ImageFileReader< InternalImageType > MovingImageReaderType;

  // read images
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[FIXED_IMG ] );
  movingImageReader->SetFileName( argv[MOVING_IMG] );

  fixedImageReader->Update(  );
  movingImageReader->Update( );

  typedef itk::TranslationTransform< double, Dimension > TransformType;
  // Software Guide : EndCodeSnippet
  typedef itk::ExhaustiveOptimizer                       OptimizerType;
  typedef itk::MattesMutualInformationImageToImageMetric<
                                    InternalImageType,
                                    InternalImageType >   MetricType;
  typedef itk::LinearInterpolateImageFunction<
                                    InternalImageType,
                                    double             > InterpolatorType;
  typedef itk::ImageRegistrationMethod< InternalImageType, InternalImageType >       RegistrationType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  TransformType::Pointer      transform     = TransformType::New();

  // Index observer (enables us to check if all positions were indeed visisted):
  IndexObserver::Pointer idxObserver = IndexObserver::New ();
  optimizer->AddObserver ( itk::IterationEvent (), idxObserver );

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  registration->SetTransform(     transform     );
  registration->SetFixedImage(    fixedImageReader->GetOutput()    );
  registration->SetMovingImage(   movingImageReader->GetOutput()   );
  registration->SetFixedImageRegion(
     fixedImageReader->GetOutput()->GetBufferedRegion() );

  typedef MetricType::ParametersType    ParametersType;


  const unsigned int spaceDimension =
                      metric->GetNumberOfParameters();

  // We start not so far from  | 2 -2 |
  ParametersType  initialPosition( spaceDimension );
  InternalImageType::PointType     origin = movingImageReader->GetOutput()->GetOrigin();
  initialPosition[0] = origin[0] ;
  initialPosition[1] = origin[1] ;
  initialPosition[2] = origin[2] ;

  optimizer->SetInitialPosition( initialPosition );


  typedef  OptimizerType::ScalesType            ScalesType;
  ScalesType    parametersScale( spaceDimension );
  InternalImageType::SpacingType   spacing = fixedImageReader->GetOutput()->GetSpacing();
  parametersScale[0] = spacing[0] ;
  parametersScale[1] = spacing[1] ;
  parametersScale[2] = spacing[2] ;

  optimizer->SetScales( parametersScale );


  optimizer->SetStepLength( 1.0 );


  typedef OptimizerType::StepsType  StepsType;
  StepsType steps( Dimension );
  steps[0] = 10;
  steps[1] = 10;
  steps[1] = 10;

  optimizer->SetNumberOfSteps( steps );

  registration->SetInitialTransformParameters( transform->GetParameters() );

  try
    {
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error occurred during Optimization" << std::endl;
    std::cout << "Location    = " << e.GetLocation()    << std::endl;
    std::cout << "Description = " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }


  bool minimumValuePass = vnl_math_abs ( optimizer->GetMinimumMetricValue() - -10 ) < 1E-3;

  std::cout << "MinimumMetricValue = " << optimizer->GetMinimumMetricValue() << std::endl;
  std::cout << "Minimum Position = " << optimizer->GetMinimumMetricValuePosition() << std::endl;

  bool maximumValuePass = vnl_math_abs ( optimizer->GetMaximumMetricValue() - 926 ) < 1E-3;
  std::cout << "MaximumMetricValue = " << optimizer->GetMaximumMetricValue() << std::endl;
  std::cout << "Maximum Position = " << optimizer->GetMaximumMetricValuePosition() << std::endl;

  ParametersType finalPosition = optimizer->GetMinimumMetricValuePosition();
  std::cout << "Solution        = (";
  std::cout << finalPosition[0] << ",";
  std::cout << finalPosition[1] << ")" << std::endl;

  bool visitedIndicesPass = true;
  std::vector < unsigned long > visitedIndices = idxObserver->m_VisitedIndices;

  size_t requiredNumberOfSteps = ( 2 * steps [ 0 ] + 1 ) * ( 2 * steps [ 1 ] + 1 );
  if ( visitedIndices.size () != requiredNumberOfSteps )
  {
    visitedIndicesPass = false;
  }

  std::sort ( visitedIndices.begin (), visitedIndices.end () );

  for ( size_t i = 0; i < visitedIndices.size (); ++i )
    {
    if ( visitedIndices [ i ] != i )
      {
      visitedIndicesPass = false;
      std::cout << "Mismatch in visited index " << visitedIndices [ i ] << " @ " << i << std::endl;
      break;
      }
    }

  //
  // check results to see if it is within range
  //
  bool trueParamsPass = true;
  double trueParameters[2] = { 2, -2 };
  for( unsigned int j = 0; j < 2; j++ )
    {
    if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > 0.01 )
      {
      trueParamsPass = false;
      }
    }

  if( !minimumValuePass || !maximumValuePass || !trueParamsPass || !visitedIndicesPass )
    {
    std::cout << "minimumValuePass   = " << minimumValuePass << std::endl;
    std::cout << "maximumValuePass   = " << maximumValuePass << std::endl;
    std::cout << "trueParamsPass     = " << trueParamsPass << std::endl;
    std::cout << "visitedIndicesPass = " << visitedIndicesPass << std::endl;
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }


  std::cout << "Testing PrintSelf " << std::endl;
  optimizer->Print( std::cout );

  std::cout << "Test passed." << std::endl;
  return EXIT_SUCCESS;


}
