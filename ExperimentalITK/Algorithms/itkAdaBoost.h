/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkAdaBoost.h,v $
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
#ifndef __itkAdaBoost_h
#define __itkAdaBoost_h

#include "itkCSVArray2DDataObject.h"
#include "itkObject.h"
#include "itkProcessObject.h"
#include "itkVectorContainer.h"

#include <vector>

namespace itk
{

/** \class FeatureNode
 *
 */

template<class TReal = float>
class ITK_EXPORT FeatureNode
: public Object
{
public:

  /** Standard class typedefs. */
  typedef FeatureNode              Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef TReal                    RealType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( FeatureNode, Object );

  /** Define the sign type designating binary class membership */
  enum MembershipSignType { BACKGROUND = -1, FOREGROUND = +1 };

  /** Set/Get the feature value **/
  itkSetMacro( Value, RealType );
  itkGetConstMacro( Value, RealType );

  /** Set/Get the membership sign (BACKGROUND=-1, FOREGROUND=1) */
  itkSetMacro( MembershipSign, MembershipSignType );
  itkGetConstMacro( MembershipSign, MembershipSignType );

  /** Set/Get the weighted rate */
  itkSetMacro( Weight, RealType );
  itkGetConstMacro( Weight, RealType );

protected:

  FeatureNode();
  virtual ~FeatureNode() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:

  RealType            m_Value;
  MembershipSignType  m_MembershipSign;
  RealType            m_Weight;
};

/** \class WeakClassifier
 *
 */

template<class TFeatureNode>
class ITK_EXPORT WeakClassifier
: public Object
{
public:

  /** Standard class typedefs. */
  typedef WeakClassifier           Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( WeakClassifier, Object );

  typedef TFeatureNode                                  FeatureNodeType;
  typedef typename FeatureNodeType::Pointer             FeatureNodePointer;
  typedef typename FeatureNodeType::RealType            RealType;
  typedef typename FeatureNodeType::MembershipSignType  MembershipSignType;
  typedef unsigned long                                 SizeType;
  typedef std::vector<FeatureNodePointer>               SingleFeatureObservationsContainerType;

  /** Set/Get the feature ID **/
  itkSetMacro( FeatureID, SizeType );
  itkGetConstMacro( FeatureID, SizeType );

  /** Set/Get the membership sign (BACKGROUND=-1, FOREGROUND=1) */
  itkSetMacro( MembershipSign, MembershipSignType );
  itkGetConstMacro( MembershipSign, MembershipSignType );

  /** Set/Get the threshold */
  itkSetMacro( Threshold, RealType );
  itkGetConstMacro( Threshold, RealType );

  /** Get the weighted rate */
  itkSetMacro( WeightedRate, RealType );
  itkGetConstMacro( WeightedRate, RealType );

  /** Set/Get the true error */
  itkSetMacro( TrueError, RealType );
  itkGetConstMacro( TrueError, RealType );

  /** Set/Get the weighted error */
  itkSetMacro( WeightedError, RealType );
  itkGetConstMacro( WeightedError, RealType );

  /** Set/Get the observations for a single feature */
  itkSetMacro( SingleFeatureObservations, SingleFeatureObservationsContainerType );
  itkGetConstMacro( SingleFeatureObservations, SingleFeatureObservationsContainerType );

  /** Given the single feature observations, learn the best hyperplane */
  void DoWeakLearn();

protected:

  WeakClassifier();
  virtual ~WeakClassifier() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

private:

  SizeType                                       m_FeatureID;
  MembershipSignType                             m_MembershipSign;
  RealType                                       m_Threshold;
  RealType                                       m_WeightedRate;
  RealType                                       m_TrueError;
  RealType                                       m_WeightedError;

  SingleFeatureObservationsContainerType         m_SingleFeatureObservations;
};

/** \class StrongClassifier
 *
 */

template<class TWeakClassifier>
class ITK_EXPORT StrongClassifier
: public Object
{
public:

  /** Standard class typedefs. */
  typedef StrongClassifier         Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( StrongClassifier, Object );

  /** Define the data types. */
  typedef TWeakClassifier                                        WeakClassifierType;
  typedef typename WeakClassifierType::RealType                  RealType;
  typedef typename WeakClassifierType::FeatureNodeType           FeatureNodeType;
  typedef typename WeakClassifierType::MembershipSignType        MembershipSignType;
  typedef typename WeakClassifierType::Pointer                   WeakClassifierPointer;
  typedef VectorContainer<unsigned int, WeakClassifierPointer>   ClassifierType;

  typedef std::vector<RealType>                            SingleObservationContainerType;

  typedef CSVArray2DDataObject<RealType>                   CSVObjectType;
  typedef typename CSVObjectType::Pointer                  CSVObjectPointer;
  typedef RealType                                         CSVDataType;

  /** Add a weak classifier to the solution. */
  void AddWeakClassifier( WeakClassifierType * );

  /** Clears the solution of all weak classifiers */
  void ClearWeakClassifiers();

  /** Returns membership value (and continuous membership value) */
  MembershipSignType Classify( const SingleObservationContainerType &, RealType & );

  /**
   * An input csv object is expected with the following column headers
   * \li ITERATION
   * \li FEATURE_ID
   * \li MEMBERSHIP_SIGN
   * \li THRESHOLD,
   * \li WEIGHTED_RATE
   * \li TRUE_ERROR
   * \li WEIGHTED_ERROR
   *
   * Note that the last two values are unnecessary for classification but that
   * is what is stored in the weak classifiers and saved to the csv object if
   * somebody calls \c GetClassifierCSVObject()
   */
  void GenerateClassifierFromCSVObject( const CSVObjectType * );

  /**
   * Returns a csv object with the following column headers
   * \li ITERATION
   * \li FEATURE_ID
   * \li MEMBERSHIP_SIGN
   * \li THRESHOLD,
   * \li WEIGHTED_RATE
   * \li TRUE_ERROR
   * \li WEIGHTED_ERROR
   */
  typename CSVObjectType::Pointer GetClassifierCSVObject();

protected:

  StrongClassifier();
  virtual ~StrongClassifier() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

private:

  typename ClassifierType::Pointer                 m_Classifier;
};


/** \class AdaBoost
 *
 */

template<class TStrongClassifier>
class ITK_EXPORT AdaBoost
: public Object
{
public:
  /** Standard class typedefs. */
  typedef AdaBoost                 Self;
  typedef ProcessObject            Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AdaBoost, ProcessObject );

  /** Define the data types. */
  typedef TStrongClassifier                                   StrongClassifierType;
  typedef typename StrongClassifierType::WeakClassifierType   WeakClassifierType;
  typedef typename WeakClassifierType::RealType               RealType;
  typedef typename WeakClassifierType::FeatureNodeType        FeatureNodeType;
  typedef typename WeakClassifierType::MembershipSignType     MembershipSignType;

  typedef typename WeakClassifierType::SingleFeatureObservationsContainerType SingleFeatureObservationsContainerType;
  typedef typename StrongClassifierType::SingleObservationContainerType       SingleObservationContainerType;

  /** Operations required for sorting */
  bool compare( const FeatureNodeType & node1, const FeatureNodeType & node2 )
    {
    return( node1.GetValue() < node2.GetValue() );
    }

  /** Set/Get number of iterations. */
  itkSetMacro( NumberOfIterations, unsigned int );
  itkGetConstMacro( NumberOfIterations, unsigned int );

  /** Add a single training observation (presumably with multiple features) */
  void AddTrainingObservation( MembershipSignType, SingleObservationContainerType & );

  /** Clear the set of training observations. */
  void ClearTrainingObservations();

  /** Set/Get the strong classifier. */
  itkGetConstObjectMacro( StrongClassifier, StrongClassifierType );

  void PerformTraining();

protected:
  AdaBoost();
  virtual ~AdaBoost() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:

  AdaBoost( const Self & ); // purposely not implemented
  void operator=( const Self & );          // purposely not implemented

  std::vector<SingleObservationContainerType>    m_TrainingObservations;
  std::vector<MembershipSignType>                m_MembershipSigns;

  unsigned int                                   m_NumberOfIterations;

  typename StrongClassifierType::Pointer         m_StrongClassifier;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAdaBoost.hxx"
#endif

#endif
