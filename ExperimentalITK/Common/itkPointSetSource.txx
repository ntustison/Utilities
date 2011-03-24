/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPointSetSource.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:04 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkPointSetSource_txx
#define _itkPointSetSource_txx

#include "itkPointSetSource.h"

namespace itk
{

/**
 *
 */
template<class TOutputPointSet>
PointSetSource<TOutputPointSet>
::PointSetSource()
{
  // Create the output. We use static_cast<> here because we know the default
  // output must be of type TOutputPointSet
  OutputPointSetPointer output
    = static_cast<TOutputPointSet*>(this->MakeOutput(0).GetPointer()); 

  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

  m_GenerateDataRegion = 0;
  m_GenerateDataNumberOfRegions = 0;
}

/**
 *
 */
template<class TOutputPointSet>
typename PointSetSource<TOutputPointSet>::DataObjectPointer
PointSetSource<TOutputPointSet>
::MakeOutput(unsigned int)
{
  return static_cast<DataObject*>(TOutputPointSet::New().GetPointer());
}

/**
 *
 */
template<class TOutputPointSet>
typename PointSetSource<TOutputPointSet>::OutputPointSetType *
PointSetSource<TOutputPointSet>
::GetOutput(void)
{
  if (this->GetNumberOfOutputs() < 1)
    {
    return 0;
    }
  
  return static_cast<TOutputPointSet*>
    (this->ProcessObject::GetOutput(0));
}

  
/**
 *
 */
template<class TOutputPointSet>
typename PointSetSource<TOutputPointSet>::OutputPointSetType *
PointSetSource<TOutputPointSet>
::GetOutput(unsigned int idx)
{
  return static_cast<TOutputPointSet*>
    (this->ProcessObject::GetOutput(idx));
}


/**
 *
 */
template<class TOutputPointSet>
void 
PointSetSource<TOutputPointSet>
::SetOutput(OutputPointSetType *output)
{
  itkWarningMacro(<< "SetOutput(): This method is slated to be removed from ITK.  Please use GraftOutput() in possible combination with DisconnectPipeline() instead." );
  this->ProcessObject::SetNthOutput(0, output);
}


/**
 *
 */
template<class TOutputPointSet>
void 
PointSetSource<TOutputPointSet>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
}

/**
 * 
 */
template<class TOutputPointSet>
void
PointSetSource<TOutputPointSet>
::GraftOutput(OutputPointSetType *graft)
{
  this->GraftNthOutput(0, graft);
}


/**
 * 
 */
template<class TOutputPointSet>
void
PointSetSource<TOutputPointSet>
::GraftNthOutput(unsigned int idx, DataObject *graft)
{
  if ( idx >= this->GetNumberOfOutputs() )
    {
    itkExceptionMacro(<<"Requested to graft output " << idx << 
        " but this filter only has " << this->GetNumberOfOutputs() << " Outputs.");
    }  

  if ( !graft )
    {
    itkExceptionMacro(<<"Requested to graft output that is a NULL pointer" );
    }

  DataObject * output = this->GetOutput(idx);

  // Call GraftImage to copy meta-information, regions, and the pixel container
  output->Graft( graft );
}


/**
 *
 */
template<class TOutputPointSet>
void 
PointSetSource<TOutputPointSet>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk

#endif
