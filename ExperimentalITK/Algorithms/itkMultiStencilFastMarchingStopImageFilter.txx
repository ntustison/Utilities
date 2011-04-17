/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkMultiStencilFastMarchingStopImageFilter.txx,v $
 Language:  C++
 Date:      $Date: 2010/03/07 01:50:24 $
 Version:   $Revision: 1.33 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkMultiStencilFastMarchingStopImageFilter_txx
#define __itkMultiStencilFastMarchingStopImageFilter_txx

#include "itkMultiStencilFastMarchingStopImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include <algorithm>

namespace itk {

template<class TLevelSet, class TSpeedImage>
MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::MultiStencilFastMarchingStopImageFilter() :
    m_TrialHeap() {

    this->ProcessObject::SetNumberOfRequiredInputs(0);

    OutputSizeType outputSize;
    outputSize.Fill(16);
    typename LevelSetImageType::IndexType outputIndex;
    outputIndex.Fill(0);

    m_OutputRegion.SetSize(outputSize);
    m_OutputRegion.SetIndex(outputIndex);

    m_OutputOrigin.Fill(0.0);
    m_OutputSpacing.Fill(1.0);
    m_OutputDirection.SetIdentity();
    m_OverrideOutputInformation = false;

    m_AlivePoints = NULL;
    m_TrialPoints = NULL;
    m_ProcessedPoints = NULL;

    m_SpeedConstant = 1.0;
    m_InverseSpeed = -1.0;
    m_LabelImage = LabelImageType::New();

    m_LargeValue = static_cast<PixelType> (NumericTraits<PixelType>::max()
            / 2.0);
    m_StoppingValue = static_cast<double> (m_LargeValue);
    m_CollectPoints = false;

    m_NormalizationFactor = 1.0;

    // get the index  from itksnap
    DBGINDEX[0]=188;
    DBGINDEX[1]=268;
    DBGINDEX[2]=128;

    m_TrackTimeSource = false;
}

template<class TLevelSet, class TSpeedImage>
void MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::PrintSelf(
        std::ostream& os, Indent indent) const {
    Superclass::PrintSelf(os, indent);
    os << indent << "Alive points: " << m_AlivePoints.GetPointer() << std::endl;
    os << indent << "Trial points: " << m_TrialPoints.GetPointer() << std::endl;
    os << indent << "Speed constant: " << m_SpeedConstant << std::endl;
    os << indent << "Stopping value: " << m_StoppingValue << std::endl;
    os << indent << "Large Value: " << static_cast<typename NumericTraits<
            PixelType>::PrintType> (m_LargeValue) << std::endl;
    os << indent << "Normalization Factor: " << m_NormalizationFactor
            << std::endl;
    os << indent << "Collect points: " << m_CollectPoints << std::endl;
    os << indent << "OverrideOutputInformation: ";
    os << m_OverrideOutputInformation << std::endl;
    os << indent << "OutputRegion: " << m_OutputRegion << std::endl;
    os << indent << "OutputOrigin:  " << m_OutputOrigin << std::endl;
    os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
    os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
}

template<class TLevelSet, class TSpeedImage>
void MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::GenerateOutputInformation() {

    // copy output information from input image
    Superclass::GenerateOutputInformation();

    // use user-specified output information
    if (this->GetInput() == NULL || m_OverrideOutputInformation) {
        LevelSetPointer output = this->GetOutput();
        output->SetLargestPossibleRegion(m_OutputRegion);
        output->SetOrigin(m_OutputOrigin);
        output->SetSpacing(m_OutputSpacing);
        output->SetDirection(m_OutputDirection);
    }
}

template<class TLevelSet, class TSpeedImage>
void MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::EnlargeOutputRequestedRegion(
        DataObject *output) {
    // enlarge the requested region of the output
    // to the whole data set
    TLevelSet * imgData;

    imgData = dynamic_cast<TLevelSet*> (output);
    if (imgData) {
        imgData->SetRequestedRegionToLargestPossibleRegion();
    } else {
        // Pointer could not be cast to TLevelSet *
        itkWarningMacro(<< "itk::MultiStencilFastMarchingStopImageFilter" <<
                "::EnlargeOutputRequestedRegion cannot cast "
                << typeid(output).name() << " to "
                << typeid(TLevelSet*).name() );
    }

}

template<class TLevelSet, class TSpeedImage>
void MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::Initialize(
        LevelSetImageType * output) {

    // allocate memory for the output buffer
    output->SetBufferedRegion(output->GetRequestedRegion());
    output->Allocate();

    // cache some buffered region information
    m_BufferedRegion = output->GetBufferedRegion();
    m_StartIndex = m_BufferedRegion.GetIndex();
    m_LastIndex = m_StartIndex + m_BufferedRegion.GetSize();
    typename LevelSetImageType::OffsetType offset;
    offset.Fill(1);
    m_LastIndex -= offset;

    // allocate memory for the PointTypeImage
    m_LabelImage->CopyInformation(output);
    m_LabelImage->SetBufferedRegion(output->GetBufferedRegion());
    m_LabelImage->Allocate();

    // set all output value to infinity
    typedef ImageRegionIterator<LevelSetImageType> OutputIterator;

    OutputIterator outIt(output, output->GetBufferedRegion());

    PixelType outputPixel;
    outputPixel = m_LargeValue;

    for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt) {
        outIt.Set(outputPixel);
    }

    // set all points type to FarPoint
    typedef ImageRegionIterator<LabelImageType> LabelIterator;

    LabelIterator typeIt(m_LabelImage, m_LabelImage->GetBufferedRegion());

    for (typeIt.GoToBegin(); !typeIt.IsAtEnd(); ++typeIt) {
        typeIt.Set(FarPoint);
    }

    // process input alive points
    AxisNodeType node;

    if (m_AlivePoints) {
        typename NodeContainer::ConstIterator pointsIter =
                m_AlivePoints->Begin();
        typename NodeContainer::ConstIterator pointsEnd = m_AlivePoints->End();

        for (; pointsIter != pointsEnd; ++pointsIter) {

            // get node from alive points container
            node = pointsIter.Value();

            // check if node index is within the output level set
            if (!m_BufferedRegion.IsInside(node.GetIndex())) {
                continue;
            }

            // make this an alive point
            m_LabelImage->SetPixel(node.GetIndex(), AlivePoint);

            outputPixel = node.GetValue();
            output->SetPixel(node.GetIndex(), outputPixel);

        }
    }

    // make sure the heap is empty
    while (!m_TrialHeap.empty()) {
        m_TrialHeap.pop();
    }

    // process the input trial points
    if (m_TrialPoints) {
        typename NodeContainer::ConstIterator pointsIter =
                m_TrialPoints->Begin();
        typename NodeContainer::ConstIterator pointsEnd = m_TrialPoints->End();

        for (; pointsIter != pointsEnd; ++pointsIter) {

            // get node from trial points container
            node = pointsIter.Value();

            // check if node index is within the output level set
            if (!m_BufferedRegion.IsInside(node.GetIndex())) {
                continue;
            }

            // make this a trial point
            m_LabelImage->SetPixel(node.GetIndex(), TrialPoint);

            outputPixel = node.GetValue();
            output->SetPixel(node.GetIndex(), outputPixel);

            m_TrialHeap.push(node);

        }
    }

}

template<class TLevelSet, class TSpeedImage>
void MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::GenerateData() {

    LevelSetPointer output = this->GetOutput();
    SpeedImageConstPointer speedImage = this->GetInput();

    this->Initialize(output);

    if (m_CollectPoints) {
        m_ProcessedPoints = NodeContainer::New();
    }

    if (m_TrackTimeSource) AllocateTimeSourceImage();

    // process points on the heap
    AxisNodeType node;
    double currentValue;
    double oldProgress = 0;

    this->UpdateProgress(0.0); // Send first progress event

    while (!m_TrialHeap.empty()) {

        // std::cout << "in fast marching" << std::endl;

        // get the node with the smallest value
        node = m_TrialHeap.top();
        m_TrialHeap.pop();

        // does this node contain the current value ?
        currentValue = (double) output->GetPixel(node.GetIndex());

        if (node.GetValue() != currentValue) {
            continue;
        }

        // is this node already alive ?
        if (m_LabelImage->GetPixel(node.GetIndex()) != TrialPoint) {
            continue;
        }

        //songgang: check freezing conditions
        // if need to freeze,
        //  mark the arrival time to be the stopping time
        //  remove from the heap

        // std::cout << "index="<<node.GetIndex() << " pixel=" <<m_IntensityImage->GetPixel(node.GetIndex()) << std::endl;

        if (m_IntensityImage->GetPixel(node.GetIndex())
                > m_FreezeLowestIntensity) {
            output->SetPixel(node.GetIndex(),
                    static_cast<PixelType> (m_LargeValue));
        }

        if (currentValue > m_StoppingValue) {
            break;
        }

        if (m_CollectPoints) {
            m_ProcessedPoints->InsertElement(m_ProcessedPoints->Size(), node);
        }

        // set this node as alive
        m_LabelImage->SetPixel(node.GetIndex(), AlivePoint);

        // update its neighbors
        this->UpdateNeighbors(node.GetIndex(), speedImage, output);

        // Send events every certain number of points.
        const double newProgress = currentValue / m_StoppingValue;
        if (newProgress - oldProgress > 0.01) // update every 1%
        {
            this->UpdateProgress(newProgress);
            oldProgress = newProgress;
            if (this->GetAbortGenerateData()) {
                this->InvokeEvent(AbortEvent());
                this->ResetPipeline();
                ProcessAborted e(__FILE__, __LINE__);
                e.SetDescription("Process aborted.");
                e.SetLocation(ITK_LOCATION);
                throw e;
            }
        }

    }

}

template<class TLevelSet, class TSpeedImage>
void MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::UpdateNeighbors(
        const IndexType& index, const SpeedImageType * speedImage,
        LevelSetImageType * output) {
    IndexType neighIndex = index;

    const int NeighborDef3D[][3] = {
            {-1, 0, 0},
            {1, 0, 0},
            {0, -1, 0},
            {0, 1, 0},
            {0, 0, -1},
            {0, 0, 1},
            {-1, -1, 0},
            {-1, 1, 0},
            {1 -1, 0},
            {1, 1, 0},
            {-1, 0, -1},
            {-1, 0, 1},
            {1, 0,-1},
            {1, 0, 1},
            {0, -1, -1},
            {0, -1, 1},
            {0, 1, -1},
            {0, 1, 1},
            {-1,-1,-1},
            {-1,-1,1},
            {-1,1,-1},
            {-1,1,1},
            {1,-1,-1},
            {1,-1,1},
            {1,1,-1},
            {1,1,1},
    };
    const unsigned int num_NeighborDef = 26;

    // std::cout << "looking for neighbor: index=" << index << std::endl;

    for(unsigned int m = 0; m < num_NeighborDef; m++){
        bool is_neighbor_outside = false;
        for(unsigned int j = 0; j < SetDimension; j++){
            neighIndex[j] = index[j] + NeighborDef3D[m][j];
            if (neighIndex[j] < m_StartIndex[j] || neighIndex[j] > m_LastIndex[j]) {
                is_neighbor_outside = true;
                break;
            }
        }
        if (is_neighbor_outside) continue;

        if (m_LabelImage->GetPixel(neighIndex) != AlivePoint) {

            // if ((neighIndex[0]==81-1 && neighIndex[1]==160-1 && neighIndex[2]==92-1)) std::cout << "neighIndex=" << neighIndex << std::endl;

             this->UpdateValue(neighIndex, speedImage, output);
         }
    }






//    for (unsigned int j = 0; j < SetDimension; j++) {
//        // update left neighbor
//        if (index[j] > m_StartIndex[j]) {
//            neighIndex[j] = index[j] - 1;
//        }
//        if (m_LabelImage->GetPixel(neighIndex) != AlivePoint) {
//            this->UpdateValue(neighIndex, speedImage, output);
//        }
//
//        // update right neighbor
//        if (index[j] < m_LastIndex[j]) {
//            neighIndex[j] = index[j] + 1;
//        }
//        if (m_LabelImage->GetPixel(neighIndex) != AlivePoint) {
//            this->UpdateValue(neighIndex, speedImage, output);
//        }
//
//        //reset neighIndex
//        neighIndex[j] = index[j];
//
//    }

}

template<class TLevelSet, class TSpeedImage>
void MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::AllocateTimeSourceImage(){

    typename LevelSetImageType::Pointer output = this->GetOutput();

    m_TimeSourceImage = TimeSourceImageType::New();
    m_TimeSourceImage->SetSpacing(output->GetSpacing() ) ;
    m_TimeSourceImage->SetOrigin(output->GetOrigin() );
    m_TimeSourceImage->SetDirection( output-> GetDirection() );
    m_TimeSourceImage->SetLargestPossibleRegion( output->GetLargestPossibleRegion() );
    m_TimeSourceImage->SetRequestedRegion( output->GetRequestedRegion() );
    m_TimeSourceImage->SetBufferedRegion( output->GetBufferedRegion() );
    m_TimeSourceImage->Allocate();

    std::cout << "Allocate Time Source Image " << m_TimeSourceImage << std::endl;
}




template<class TLevelSet, class TSpeedImage>
double MultiStencilFastMarchingStopImageFilter<TLevelSet, TSpeedImage>::UpdateValue(
        const IndexType& index, const SpeedImageType * speedImage,
        LevelSetImageType * output) {

    IndexType neighIndex = index;
    typename TLevelSet::PixelType neighValue;
    PixelType outputPixel;
    AxisNodeType node;

    if ((index[0]==DBGINDEX[0]-1 && index[1]==DBGINDEX[1]-1 && index[2]==DBGINDEX[2]-1))
         std::cout << "update index=" << index<<std::endl;




    //songgang: check freezing conditions
    // if need to freeze,
    //  mark the arrival time to be the stopping time
    //  remove from the heap
    if (m_IntensityImage->GetPixel(index) > m_FreezeLowestIntensity) {
        output->SetPixel(index, static_cast<PixelType> (m_LargeValue));
        m_LabelImage->SetPixel(index, TrialPoint);

        if (index[0] == DBGINDEX[0] - 1 && index[1] == DBGINDEX[1] - 1 && index[2] == DBGINDEX[2] - 1) {
            std::cout << " --------------- reach point " << index << ":"
                    << output->GetPixel(index) << std::endl;
        }

        return m_LargeValue;
    }

    // using multi stencils and choose the one with smallest (and valid) arrival time
    // note: if no real roots from 2nd order equation, choose the arrival time from the neighbor with smallest
    // arrival time

    const int StencilDef[][3][2][3] = {
            { {{-1,0,0},{1,0,0}}, {{0,-1,0},{0,1,0}}, {{0,0,-1},{0,0,1}} },
            { {{-1,0,0},{1,0,0}}, {{0,-1,1},{0,1,-1}}, {{0,-1,-1},{0,1,1}} },
            { {{0,-1,0},{0,1,0}}, {{-1,0,1},{1,0,-1}}, {{-1,0,-1},{1,0,1}} },
            { {{0,0,-1},{0,0,1}}, {{-1,1,0},{1,-1,0}}, {{-1,-1,0},{1,1,0}} },
            { {{-1,0,-1},{1,0,1}}, {{1,-1,-1},{-1,1,1}}, {{-1,-1,1},{1,1,-1}} },
            { {{-1,0,1},{1,0,-1}}, {{-1,-1,-1},{1,1,1}}, {{1,-1,1},{-1,1,-1}} },
    };

    const bool  StencilAxis[][3][3] = { //[stencil][axis in stencil][valid neighbor in global axis]
            {{1,0,0}, {0,1,0}, {0,0,1}},
            {{1,0,0}, {0,1,1}, {0,1,1}},
            {{0,1,0}, {1,0,1}, {1,0,1}},
            {{0,0,1}, {1,1,0}, {1,1,0}},
            {{1,0,1}, {1,1,1}, {1,1,1}},
            {{1,0,1}, {1,1,1}, {1,1,1}}
    };

    const unsigned int num_Stencil = 6;


    // keep track of the "best" neighbor node for current node:
    unsigned int best_stencil = -1;
    int best_axis = -1;
    int best_left_or_right = -1;

    double solution_all_stencil = m_LargeValue;
    for(unsigned int s=0; s<num_Stencil; s++) { // iterate over all stencils
        for(unsigned int j=0; j<SetDimension; j++) { // iterate over the 3 axis for each stencil

            node.SetValue(m_LargeValue);

            for(unsigned int k=0; k<2; k++) {// left and right of each axis

                bool is_neighbor_outside = false;
                for(unsigned int q=0; q<SetDimension; q++){
                    neighIndex[q] = index[q] + StencilDef[s][j][k][q];
                    if (neighIndex[q] > m_LastIndex[q] || neighIndex[q] < m_StartIndex[j]) {
                        is_neighbor_outside = true; break;
                    }
                }

                if ((index[0]==DBGINDEX[0]-1 && index[1]==DBGINDEX[1]-1 && index[2]==DBGINDEX[2]-1))
                    std::cout << "index=" << index<<" neighIndex="<<neighIndex<<" neighValue="<<output->GetPixel(neighIndex) <<std::endl;


                if (is_neighbor_outside) continue;


                if (m_LabelImage->GetPixel(neighIndex) == AlivePoint) {
                       outputPixel = output->GetPixel(neighIndex);
                       neighValue = outputPixel;

                       if (node.GetValue() > neighValue) {
                           node.SetValue(neighValue);
                           node.SetIndex(neighIndex);
                           node.SetLeftOrRight(k);
                       }
                }
            }

            // put the minimum neighbor onto the heap
            m_NodesUsed[j] = node;
            m_NodesUsed[j].SetAxis(j);

        } // end of each axis in each stencil

        // construct the quadratic equation
        // sort the local list
        std::sort(m_NodesUsed, m_NodesUsed + SetDimension);

        // solve quadratic equation
        double aa, bb, cc;
        double solution = m_LargeValue;

        aa = 0.0;
        bb = 0.0;
        if (speedImage) {
            typedef typename SpeedImageType::PixelType SpeedPixelType;
            cc = (double) speedImage->GetPixel(index) / m_NormalizationFactor;
            cc = -1.0 * vnl_math_sqr(1.0 / cc);
        } else {
            cc = m_InverseSpeed;
        }

        OutputSpacingType spacing = this->GetOutput()->GetSpacing();

        double discrim;

        for (unsigned int j = 0; j < SetDimension; j++) {
            node = m_NodesUsed[j];

            if (solution >= node.GetValue()) {
                const int axis = node.GetAxis();
                // const double spaceFactor = vnl_math_sqr(1.0 / spacing[axis]);
                double spaceFactor = 0.0;
                for(unsigned int p=0; p<SetDimension; p++)
                    if (StencilAxis[s][axis][p]) spaceFactor += vnl_math_sqr(spacing[p]);
                spaceFactor = 1.0 / spaceFactor;

                const double value = double(node.GetValue());
                aa += spaceFactor;
                bb += value * spaceFactor;
                cc += vnl_math_sqr(value) * spaceFactor;


                discrim = vnl_math_sqr(bb) - aa * cc;
                if (discrim < 0.0) {
                    // Discriminant of quadratic eqn. is negative
                    ExceptionObject err(__FILE__, __LINE__);
                    err.SetLocation(ITK_LOCATION);
                    err.SetDescription(
                            "Discriminant of quadratic equation is negative");
                    // throw err;
                    continue;
                }

                double solution1 = (vcl_sqrt(discrim) + bb) / aa;
                if (solution > solution1) solution = solution1;

                if ((index[0]==DBGINDEX[0]-1 && index[1]==DBGINDEX[1]-1 && index[2]==DBGINDEX[2]-1))
                     std::cout << "aa/bb/cc" << aa << " " << bb <<" " << cc << " solution1=" << solution1 << std::endl;


            } else {
                break;
            }
        }

        if (solution_all_stencil > solution) {
            solution_all_stencil = solution;

            if (m_TrackTimeSource) {
                best_stencil = s;
                best_axis = m_NodesUsed[0].GetAxis();
                best_left_or_right = m_NodesUsed[0].GetLeftOrRight();

                // int current_time_source[3];
                TimeSourceVectorType current_time_source;
                for(unsigned int bbb = 0; bbb < 3; bbb++)
//                    current_time_source[bbb] = 70 * (StencilDef[best_stencil][best_axis][best_left_or_right][bbb] + 2); // -1,0,1 -> 70, 140, 210 (RGB value)
                    current_time_source[bbb] = StencilDef[best_stencil][best_axis][best_left_or_right][bbb]; // -1,0,1 -> 70, 140, 210 (RGB value)
//                current_time_source[0] = 71;
//                current_time_source[1] = 53;
//                current_time_source[2] = 146;

                m_TimeSourceImage->SetPixel(index, current_time_source);
            }
        }

    } // end iteration of all stencils

    if ((index[0]==DBGINDEX[0]-1 && index[1]==DBGINDEX[1]-1 && index[2]==DBGINDEX[2]-1))
         std::cout << "solution_all_stencil=" << solution_all_stencil << std::endl;


    if (solution_all_stencil < m_LargeValue) {
        // write solution to m_OutputLevelSet
        outputPixel = static_cast<PixelType> (solution_all_stencil);
        output->SetPixel(index, outputPixel);

        // insert point into trial heap
        m_LabelImage->SetPixel(index, TrialPoint);
        node.SetValue(static_cast<PixelType> (solution_all_stencil));
        node.SetIndex(index);
        m_TrialHeap.push(node);

        if (index[0] == DBGINDEX[0] - 1 && index[1] == DBGINDEX[1] - 1 && index[2] == DBGINDEX[2] - 1) {
            std::cout << " --------------- why here point??? " << index << ":"
                    << output->GetPixel(index) << " i="
                    << m_IntensityImage->GetPixel(index) << " t="
                    << m_FreezeLowestIntensity << std::endl;
        }

    }

    return solution_all_stencil;







//    // sort the local list
//    std::sort(m_NodesUsed, m_NodesUsed + SetDimension);
//
//    // solve quadratic equation
//    double aa, bb, cc;
//    double solution = m_LargeValue;
//
//    aa = 0.0;
//    bb = 0.0;
//    if (speedImage) {
//        typedef typename SpeedImageType::PixelType SpeedPixelType;
//        cc = (double) speedImage->GetPixel(index) / m_NormalizationFactor;
//        cc = -1.0 * vnl_math_sqr(1.0 / cc);
//    } else {
//        cc = m_InverseSpeed;
//    }
//
//    OutputSpacingType spacing = this->GetOutput()->GetSpacing();
//
//    double discrim;
//
//    for (unsigned int j = 0; j < SetDimension; j++) {
//        node = m_NodesUsed[j];
//
//        if (solution >= node.GetValue()) {
//            const int axis = node.GetAxis();
//            const double spaceFactor = vnl_math_sqr(1.0 / spacing[axis]);
//            const double value = double(node.GetValue());
//            aa += spaceFactor;
//            bb += value * spaceFactor;
//            cc += vnl_math_sqr(value) * spaceFactor;
//
//            discrim = vnl_math_sqr(bb) - aa * cc;
//            if (discrim < 0.0) {
//                // Discriminant of quadratic eqn. is negative
//                ExceptionObject err(__FILE__, __LINE__);
//                err.SetLocation(ITK_LOCATION);
//                err.SetDescription(
//                        "Discriminant of quadratic equation is negative");
//                throw err;
//            }
//
//            solution = (vcl_sqrt(discrim) + bb) / aa;
//        } else {
//            break;
//        }
//    }
//
//    if (solution < m_LargeValue) {
//        // write solution to m_OutputLevelSet
//        outputPixel = static_cast<PixelType> (solution);
//        output->SetPixel(index, outputPixel);
//
//        // insert point into trial heap
//        m_LabelImage->SetPixel(index, TrialPoint);
//        node.SetValue(static_cast<PixelType> (solution));
//        node.SetIndex(index);
//        m_TrialHeap.push(node);
//
//        if (index[0] == 62 - 1 && index[1] == 12 - 1 && index[2] == 94 - 1) {
//            std::cout << " --------------- why here point??? " << index << ":"
//                    << output->GetPixel(index) << " i="
//                    << m_IntensityImage->GetPixel(index) << " t="
//                    << m_FreezeLowestIntensity << std::endl;
//        }
//
//    }
//
//    return solution;

}

} // namespace itk


#endif
