#include "itkAdaBoost.h"

int main( int argc, char *argv[] )
{
  typedef float                                        RealType;
  typedef itk::FeatureNode<RealType>                   FeatureNodeType;
  typedef FeatureNodeType::MembershipSignType          MembershipSignType;
  typedef itk::WeakClassifier<FeatureNodeType>         WeakClassifierType;
  typedef itk::StrongClassifier<WeakClassifierType>    StrongClassifierType;

  typedef itk::AdaBoost<StrongClassifierType>           AdaBoostType;
  typedef AdaBoostType::SingleObservationContainerType ObservationContainerType;


  AdaBoostType::Pointer adaboost = AdaBoostType::New();
  adaboost->SetNumberOfIterations( 3 );

  RealType xvalues[] = { 0, 1, 2, 3, 4, 5 };
  MembershipSignType yvalues[] = { FeatureNodeType::FOREGROUND,
                                   FeatureNodeType::FOREGROUND,
                                   FeatureNodeType::FOREGROUND,
                                   FeatureNodeType::BACKGROUND,
                                   FeatureNodeType::BACKGROUND,
                                   FeatureNodeType::BACKGROUND };

  for( unsigned int d = 0; d < 10; d++ )
    {
    ObservationContainerType observation;
    observation.push_back( xvalues[d] );

    adaboost->AddTrainingObservation( yvalues[d], observation );
    }

  adaboost->PerformTraining();

  adaboost->GetStrongClassifier()->Print( std::cout, 3 );

  return EXIT_SUCCESS;
}
