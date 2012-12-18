/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkAdaBoost.hxx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAdaBoost_hxx
#define __itkAdaBoost_hxx

#include "itkAdaBoost.h"

#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"

#include <algorithm>

namespace itk
{

/**********************************************************/
/*               FeatureNode class definitions            */
/**********************************************************/

template<class TReal>
FeatureNode<TReal>
::FeatureNode() :
  m_Value( 0.0 ),
  m_MembershipSign( FOREGROUND  ),
  m_Weight( 0.0 )
{}

template<class TReal>
void
FeatureNode<TReal>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  os << indent << "Value:      " << this->m_Value << std::endl;
  os << indent << "Membership Sign: " << this->m_MembershipSign << std::endl;
  os << indent << "Weight:      " << this->m_Weight << std::endl;
}

/**********************************************************/
/*               WeakClassifier class definitions         */
/**********************************************************/

template<class TFeatureNode>
WeakClassifier<TFeatureNode>
::WeakClassifier() :
  m_FeatureID( 0 ),
  m_MembershipSign( FeatureNodeType::FOREGROUND ),
  m_Threshold( 0.0 ),
  m_WeightedRate( 0.0 ),
  m_TrueError( 0.0 ),
  m_WeightedError( 0.0 )
{
  this->m_SingleFeatureObservations.clear();
}

template<class TFeatureNode>
void
WeakClassifier<TFeatureNode>
::DoWeakLearn()
{
  /** WeakLearn */

  if( this->m_SingleFeatureObservations.size() == 0 )
    {
    itkExceptionMacro( "There are no observations for this feature." );
    }

  RealType sumOfWeights = 0.0;
  RealType sumOfForegroundWeights = 0.0;

  typename SingleFeatureObservationsContainerType::const_iterator it =
    this->m_SingleFeatureObservations.begin();
  while( it != this->m_SingleFeatureObservations.end() )
    {
    if( ( *it )->GetMembershipSign() == FeatureNodeType::FOREGROUND )
      {
      sumOfForegroundWeights += ( *it )->GetWeight();
      }
    ++it;
    }

  RealType sumOfBackgroundWeights = 1.0 - sumOfForegroundWeights;

  this->m_Threshold = ( this->m_SingleFeatureObservations[0] )->GetValue();
  if( sumOfForegroundWeights >= sumOfBackgroundWeights )
    {
    this->m_MembershipSign = FeatureNodeType::FOREGROUND;
    sumOfWeights = sumOfForegroundWeights;
    }
  else
    {
    this->m_MembershipSign = FeatureNodeType::BACKGROUND;
    sumOfWeights = sumOfBackgroundWeights;
    }

  it = this->m_SingleFeatureObservations.begin();
  while( it != this->m_SingleFeatureObservations.end() )
    {
    typename FeatureNodeType::Pointer currentNode = *it;
    unsigned int index = it - this->m_SingleFeatureObservations.begin();
    if( index < this->m_SingleFeatureObservations.size() - 1 )
      {
      typename FeatureNodeType::Pointer nextNode = this->m_SingleFeatureObservations[index + 1];
      while( currentNode->GetValue() == nextNode->GetValue() )
        {
        if( currentNode->GetMembershipSign() == FeatureNodeType::FOREGROUND )
          {
          sumOfForegroundWeights -= currentNode->GetWeight();
          }
        else
          {
          sumOfForegroundWeights += currentNode->GetWeight();
          }
        ++it;
        currentNode = *it;
        index = it - this->m_SingleFeatureObservations.begin();
        nextNode = this->m_SingleFeatureObservations[index + 1];
        }
      }
    if( currentNode->GetMembershipSign() == FeatureNodeType::FOREGROUND )
      {
      sumOfForegroundWeights -= currentNode->GetWeight();
      }
    else
      {
      sumOfForegroundWeights += currentNode->GetWeight();
      }

    sumOfBackgroundWeights = 1.0 - sumOfForegroundWeights;

    if( sumOfForegroundWeights > sumOfWeights )
      {
      sumOfWeights = sumOfForegroundWeights;
      this->m_Threshold = currentNode->GetValue();
      this->m_MembershipSign = FeatureNodeType::FOREGROUND;
      }
    if( sumOfBackgroundWeights > sumOfWeights )
      {
      sumOfWeights = sumOfBackgroundWeights;
      this->m_Threshold = currentNode->GetValue();
      this->m_MembershipSign = FeatureNodeType::BACKGROUND;
      }
    ++it;
    }

  this->m_WeightedRate = sumOfWeights;

  std::cout << this->m_Threshold << " -> " << this->m_WeightedRate << std::endl;
}

template<class TFeatureNode>
void
WeakClassifier<TFeatureNode>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  os << indent << "Feature ID:              " << this->m_FeatureID << std::endl;
  os << indent << "Membership Sign:         " << this->m_MembershipSign << std::endl;
  os << indent << "Threshold:               " << this->m_Threshold << std::endl;
  os << indent << "Weighted Rate:           " << this->m_WeightedRate << std::endl;
  os << indent << "True Error:              " << this->m_TrueError << std::endl;
  os << indent << "Weighted Error:          " << this->m_WeightedError << std::endl;
  os << indent << "Number of observations:  " << this->m_SingleFeatureObservations.size() << std::endl;
}

/**********************************************************/
/*             StrongClassifier class definitions         */
/**********************************************************/

template<class TWeakClassifier>
StrongClassifier<TWeakClassifier>
::StrongClassifier()
{
  this->m_Classifier = ClassifierType::New();
  this->m_Classifier->Initialize();
}

template<class TWeakClassifier>
void
StrongClassifier<TWeakClassifier>
::AddWeakClassifier( WeakClassifierType * weakClassifier )
{
  typename ClassifierType::ElementIdentifier index = this->m_Classifier->Size();
  this->m_Classifier->InsertElement( index, weakClassifier );
  this->Modified();
}

template<class TWeakClassifier>
void
StrongClassifier<TWeakClassifier>
::ClearWeakClassifiers()
{
  this->m_Classifier->Initialize();
  this->Modified();
}

template<class TWeakClassifier>
typename StrongClassifier<TWeakClassifier>::MembershipSignType
StrongClassifier<TWeakClassifier>
::Classify( const SingleObservationContainerType & observation, RealType & continuousHypothesis )
{
  if( this->m_Classifier->Size() == 0 )
    {
    itkExceptionMacro( "The number of weak classifiers is 0." );
    }

  if( this->m_Classifier->Size() > observation.size() )
    {
    itkExceptionMacro( "The number of weak classifiers exceeds the number of features in the observation."
      << "This indicates a mismatch between the training and specification of the individual features." );
    }

  MembershipSignType membershipSign = FeatureNodeType::FOREGROUND;
  continuousHypothesis = 0.0;

  typename ClassifierType::ConstIterator It;
  for( It = this->m_Classifier()->Begin(); It != this->m_Classifier()->End(); ++It )
    {
    typename WeakClassifierType::Pointer weakClassifier = It.Value();

    RealType preFactor = -1.0;
    if( ( weakClassifier->GetMembershipSign() == FeatureNodeType::BACKGROUND &&
      observation[It.Index()] <= weakClassifier->GetThreshold() ) ||
      ( weakClassifier->GetMembershipSign() == FeatureNodeType::FOREGROUND &&
      observation[It.Index()] > weakClassifier->GetThreshold() ) )
      {
      preFactor = 1.0;
      }

    RealType alpha = 0.5 * vcl_log( weakClassifier->GetWeightedRate() / ( 1.0 - weakClassifier->GetWeightedRate() ) );
    continuousHypothesis += ( preFactor * alpha );
    }

  if( continuousHypothesis > 0.5 )
    {
    membershipSign = FeatureNodeType::FOREGROUND;
    }
  else
    {
    membershipSign = FeatureNodeType::BACKGROUND;
    }

  return membershipSign;
}

template<class TWeakClassifier>
void
StrongClassifier<TWeakClassifier>
::GenerateClassifierFromCSVObject( const CSVObjectType * csvClassifier )
{
  // We assume a csv file with the following column headers:
  // ITERATION,FEATURE_ID,MEMBERSHIP_SIGN,THRESHOLD,WEIGHTED_RATE,TRUE_ERROR,WEIGHTED_ERROR

  if( ( csvClassifier->GetMatrix() ).cols() != 7 )
    {
    itkExceptionMacro( "CSV object assumes 7 columns: ITERATION,FEATURE_ID,MEMBERSHIP_SIGN,THRESHOLD,WEIGHTED_RATE,TRUE_ERROR,WEIGHTED_ERROR." );
    }

  this->m_Classifiers->Initialize();

  for( unsigned int i = 0; i < ( csvClassifier->GetMatrix() ).rows(); ++i )
    {
    typename WeakClassifierType::Pointer weakClassifier = WeakClassifierType::New();
    weakClassifier->SetFeatureID( csvClassifier->GetMatrixData( i, 1 ) );
    weakClassifier->SetMembershipSign( csvClassifier->GetMatrixData( i, 2 ) );
    weakClassifier->SetThreshold( csvClassifier->GetMatrixData( i, 3 ) );
    weakClassifier->SetWeightedRate( csvClassifier->GetMatrixData( i, 4 ) );
    weakClassifier->SetTrueError( csvClassifier->GetMatrixData( i, 5 ) );
    weakClassifier->SetWeightedError( csvClassifier->GetMatrixData( i, 6 ) );

    unsigned int iteration = static_cast<unsigned int>( csvClassifier->GetMatrixData( i, 0 ) );

    this->m_Classifier->InsertElement( iteration, weakClassifier );
    }

  this->m_Classifier->Squeeze();
}

template<class TWeakClassifier>
typename StrongClassifier<TWeakClassifier>::CSVObjectPointer
StrongClassifier<TWeakClassifier>
::GetClassifierCSVObject()
{
  // We assume a csv file with the following column headers:
  // ITERATION,FEATURE_ID,MEMBERSHIP_SIGN,THRESHOLD,WEIGHTED_RATE,TRUE_ERROR,WEIGHTED_ERROR

  typename CSVObjectType::Pointer csvClassifier = CSVObjectType::New();

  csvClassifier->SetMatrixSize( this->m_Classifier()->Size(), 7 );

  csvClassifier->RowHeadersPushBack( "ITERATION" );
  csvClassifier->RowHeadersPushBack( "FEATURE_ID" );
  csvClassifier->RowHeadersPushBack( "MEMBERSHIP_SIGN" );
  csvClassifier->RowHeadersPushBack( "THRESHOLD" );
  csvClassifier->RowHeadersPushBack( "WEIGHTED_RATE" );
  csvClassifier->RowHeadersPushBack( "TRUE_ERROR" );
  csvClassifier->RowHeadersPushBack( "WEIGHTED_ERROR" );

  typename ClassifierType::ConstIterator It;
  for( It = this->m_Classifier()->Begin(); It != this->m_Classifier()->End(); ++It )
    {
    typename WeakClassifierType::Pointer weakClassifier = It.Value();
    unsigned int iteration = It.Index();

    csvClassifier->SetMatrixData( iteration, 0, static_cast<CSVDataType>( iteration ) );
    csvClassifier->SetMatrixData( iteration, 1, static_cast<CSVDataType>( weakClassifier->GetFeatureID() ) );
    csvClassifier->SetMatrixData( iteration, 2, static_cast<CSVDataType>( weakClassifier->GetMembershipSign() ) );
    csvClassifier->SetMatrixData( iteration, 3, static_cast<CSVDataType>( weakClassifier->GetThreshold() ) );
    csvClassifier->SetMatrixData( iteration, 4, static_cast<CSVDataType>( weakClassifier->GetWeightedRate() ) );
    csvClassifier->SetMatrixData( iteration, 5, static_cast<CSVDataType>( weakClassifier->GetTrueError() ) );
    csvClassifier->SetMatrixData( iteration, 6, static_cast<CSVDataType>( weakClassifier->GetWeightedError() ) );
    }

  return csvClassifier;
}

template<class TWeakClassifier>
void
StrongClassifier<TWeakClassifier>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  os << indent << "Number of weak classifiers comprising the strong classifier: " <<  this->m_Classifier->Size() << std::endl;
}

/**********************************************************/
/*             AdaBoost class definitions                 */
/**********************************************************/

template<class TStrongClassifier>
AdaBoost<TStrongClassifier>
::AdaBoost()
{
  this->m_TrainingObservations.clear();
  this->m_MembershipSigns.clear();

  this->m_StrongClassifier = StrongClassifierType::New();
}

template<class TStrongClassifier>
void
AdaBoost<TStrongClassifier>
::AddTrainingObservation( MembershipSignType membershipSign, SingleObservationContainerType & observation )
{
  if( this->m_TrainingObservations.size() > 0 && observation.size() != this->m_TrainingObservations[0].size() )
    {
    itkExceptionMacro( "The number of features for each observation must be the same." );
    }
  else
    {
    this->m_TrainingObservations.push_back( observation );
    this->m_MembershipSigns.push_back( membershipSign );
    this->Modified();
    }
}

template<class TStrongClassifier>
void
AdaBoost<TStrongClassifier>
::ClearTrainingObservations()
{
  if( this->m_TrainingObservations.size() > 0 )
    {
    this->m_TrainingObservations.clear();
    this->m_MembershipSigns.clear();
    this->m_Modified();
    }
}

template<class TStrongClassifier>
void
AdaBoost<TStrongClassifier>
::PerformTraining()
{
  unsigned int numberOfObservations = this->m_TrainingObservations.size();

  if( numberOfObservations == 0 )
    {
    return;
    }

  /** Check to see if we have two labels/classes/memberships in the set of observations */

  unsigned int numberOfForegroundObservations = 0;
  unsigned int numberOfBackgroundObservations = 0;

  typename std::vector<MembershipSignType>::const_iterator itM;
  for( itM = this->m_MembershipSigns.begin(); itM != this->m_MembershipSigns.end(); ++itM )
    {
    if( *itM == FeatureNodeType::FOREGROUND )
      {
      numberOfForegroundObservations++;
      }
    else
      {
      numberOfBackgroundObservations++;
      }
    }

  if( numberOfBackgroundObservations == 0 || numberOfForegroundObservations == 0 )
    {
    itkExceptionMacro( "All training observations are from a single membership." );
    }

  vnl_vector<RealType> weights( numberOfObservations, 1.0 / static_cast<RealType>( numberOfObservations ) );

  /** Train */

  unsigned int numberOfFeatures = this->m_TrainingObservations[0].size();

  for( unsigned int i = 0; i < this->m_NumberOfIterations; i++ )
    {
    typename WeakClassifierType::Pointer optimalWeakClassifier = WeakClassifierType::New();

    for( unsigned int j = 0; j < numberOfFeatures; j++ )
      {
      SingleFeatureObservationsContainerType singleFeatureObservations;

      typename std::vector<SingleObservationContainerType>::const_iterator it;
      for( it = this->m_TrainingObservations.begin(); it != this->m_TrainingObservations.end(); ++it )
        {
        unsigned int index = it - this->m_TrainingObservations.begin();

        typename FeatureNodeType::Pointer node = FeatureNodeType::New();
        node->SetValue( ( *it )[j] );
        node->SetMembershipSign( this->m_MembershipSigns[index] );
        node->SetWeight( weights[index] );

        singleFeatureObservations.push_back( node );
        }
      std::sort( singleFeatureObservations.begin(), singleFeatureObservations.end() );

      typename WeakClassifierType::Pointer weakClassifier = WeakClassifierType::New();
      weakClassifier->SetFeatureID( j );
      weakClassifier->SetSingleFeatureObservations( singleFeatureObservations );
      weakClassifier->DoWeakLearn();

      if( weakClassifier->GetWeightedRate() > optimalWeakClassifier->GetWeightedRate() )
        {
        optimalWeakClassifier->SetMembershipSign( weakClassifier->GetMembershipSign() );
        optimalWeakClassifier->SetWeightedRate( weakClassifier->GetWeightedRate() );
        optimalWeakClassifier->SetThreshold( weakClassifier->GetThreshold() );
        optimalWeakClassifier->SetFeatureID( j );
        }
      }

    RealType optimalWeightedRate = optimalWeakClassifier->GetWeightedRate();
    MembershipSignType optimalMembershipSign = optimalWeakClassifier->GetMembershipSign();
    RealType optimalThreshold = optimalWeakClassifier->GetThreshold();

    RealType alpha = 0.5 * vcl_log( optimalWeightedRate / ( 1.0 - optimalWeightedRate ) );
    RealType weightedError = 1.0;
    RealType trueError = 0.0;
    RealType CH = 0.0;              // real value either -1 or +1
    RealType H = 0.0;
    RealType F = 0.0;

    typename std::vector<SingleObservationContainerType>::const_iterator it;
    for( it = this->m_TrainingObservations.begin(); it != this->m_TrainingObservations.end(); ++it )
      {
      unsigned int index = it - this->m_TrainingObservations.begin();

      if( ( optimalMembershipSign == FeatureNodeType::FOREGROUND && ( *it )[index] > optimalThreshold ) ||
        ( optimalMembershipSign == FeatureNodeType::BACKGROUND && ( *it )[index] <= optimalThreshold ) )
        {
        H = 1.0;
        }
      else
        {
        H = -1.0;
        }

      if( H == static_cast<RealType>( this->m_MembershipSigns[index] ) )
        {
        weightedError -= weights[index];
        }


      CH += alpha * H;

      if( CH * static_cast<RealType>( this->m_MembershipSigns[index] ) < 0.0 )
        {
        trueError += 1.0 / static_cast<RealType>( numberOfObservations );
        }
      weights[index] *= vcl_exp( -alpha * H * static_cast<RealType>( this->m_MembershipSigns[index] ) );
      }

    std::cout << i << ": " << optimalThreshold << ", " << weightedError << ", " << trueError << ", " << CH << std::endl;

    weights /= weights.sum();


    optimalWeakClassifier->SetWeightedError( weightedError );
    optimalWeakClassifier->SetTrueError( trueError );

    this->m_StrongClassifier->AddWeakClassifier( optimalWeakClassifier );
    }
}

template<class TStrongClassifier>
void
AdaBoost<TStrongClassifier>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  os << indent << "Number of iterations:               " << this->m_NumberOfIterations << std::endl;
  os << indent << "Number of observations:             " << this->m_TrainingObservations.size() << std::endl;

  if( this->m_TrainingObservations.size() > 0 )
    {
    os << indent << "Number of features per observation: " << this->m_TrainingObservations[0].size() << std::endl;
    }
}

} // end namespace itk
#endif
