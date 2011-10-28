// BSpline testng program
#include <iostream>
#include <itkImage.h>
//#include <itkBSplineDeformableTransformInitializer.h>
#include <itkBSplineDeformableTransform.h>
//#include <itkBSplineTransformInitializer.h>
#include <itkBSplineTransform.h>
#include <itkLinearInterpolateImageFunction.h>

const unsigned int dim=2;
const unsigned int order=3;
// params [6, 6, 0.625, -5.70833, 1.75, 2.08333, 1, 0, 0, 1]
typedef itk::Image<unsigned char, dim> ImageType;
typedef itk::BSplineDeformableTransform<double, dim, order> TransformType1;
//typedef itk::BSplineDeformableTransformInitializer<TransformType,ImageType> TransformInitializerType1;
typedef itk::BSplineTransform<double, dim, order> TransformType2;
//typedef itk::BSplineTransformInitializer<TransformType,ImageType> TransformInitializerType2;


template<class ImageType>
typename ImageType::PointType GetImageCorner(typename ImageType::Pointer image, unsigned int c, double step=0) {
   const unsigned int dim=ImageType::ImageDimension;
   typename ImageType::SizeType sz=image->GetLargestPossibleRegion().GetSize();
   typename itk::ContinuousIndex<double,2> cornerIndex;
   int pow2=1;
   for(unsigned int d=0;d<dim;++d) {
      cornerIndex[d]=((c & pow2)>0)?(static_cast<double>(sz[d])-1.0-step):0.0+step;
      pow2*=2;
   }
   typename ImageType::PointType p;
   image->TransformContinuousIndexToPhysicalPoint(cornerIndex, p);
   return p;
}


int main(int argc, char** argv) {
   try {
   TransformType1::Pointer transform1=TransformType1::New();
   TransformType2::Pointer transform2=TransformType2::New();

/*   ImageType::Pointer image=ImageType::New();
   ImageType::IndexType ix;
   ImageType::SizeType sz;
   sz[0]=5; sz[1]=6; ix.Fill(0);
   ImageType::RegionType region(ix,sz);
   image->SetRegions(region);
   ImageType::PointType origin;
   origin[0]=3;origin[1]=-3;
   image->SetOrigin(origin);
   ImageType::SpacingType spacing;
   spacing.Fill(1);
   image->SetSpacing(spacing);
   ImageType::DirectionType direction;
   direction.SetIdentity();
   image->SetDirection(direction);


   TransformInitializerType::Pointer anInitializer1=TransformInitializerType1::New();

   anInitializer1->SetImage(image);
   anInitializer1->SetTransform(transform);
   TransformInitializerType::TransformType::SizeType aTransformSize;
   for(unsigned int i=0;i<dim;++i) {aTransformSize[i]=3;}//+i;}
#if ITK_VERSION_MAJOR==3
   anInitializer1->SetGridSizeInsideTheImage(aTransformSize);
#else
   transform->SetTransformDomainMeshSize(aTransformSize);
#endif
   anInitializer1->InitializeTransform();
   std::cout << "For reference, the image corners are: " << std::endl;
   for(unsigned int i=0;i<4;++i) { std::cout << "   " << i << " : " << GetImageCorner<ImageType>(image,i) << std::endl;}
   */
   std::cout << "The transform is a: " << transform1->GetNameOfClass() << std::endl;
   std::cout << "--- " << transform1->GetTransformTypeAsString() << std::endl;
   std::cout << "The transform is a: " << transform2->GetNameOfClass() << std::endl;
   std::cout << "--- " << transform2->GetTransformTypeAsString() << std::endl;

   TransformType1::ParametersType fparam=transform1->GetFixedParameters();
   fparam[0]=6;
   fparam[1]=6;
   fparam[2]=0.625;
   fparam[3]=-5.70833;
   fparam[4]=1.75;
   fparam[5]=2.08333;
   fparam[6]=1;
   fparam[7]=0;
   fparam[8]=0;
   fparam[9]=1;
   transform1->SetFixedParameters(fparam);
   transform2->SetFixedParameters(fparam);

   // wierdly, the fixed parameters dont "stick" in the BSplineTransform
   std::cout << "The desired fixed parameters are: " << fparam << std::endl;
   std::cout << "Fixed parameters 1: " << transform1->GetFixedParameters() << std::endl;
   std::cout << "Fixed parameters 2: " << transform2->GetFixedParameters() << std::endl;

   // Set The parameters to some sort of radial bump
   std::cout << "Number of parameters 1: " << transform1->GetNumberOfParameters() << std::endl;
   std::cout << "Number of parameters 2: " << transform2->GetNumberOfParameters() << std::endl;
   TransformType1::ParametersType param=transform1->GetParameters();

   for(unsigned int r=0;r<fparam[0];++r) {
      for(unsigned int c=0;c<fparam[1];++c) {
         double dx=r-fparam[0]/2;
         double dy=c-fparam[1]/2;
         double xval=5.0 / (dx<0?-1+dx:1+dx);
         double yval=4.0 / (dy<0?-1+dy:1+dy);
         param[r*fparam[1]+c]=xval;
         param[fparam[0]*fparam[1]+r*fparam[1]+c]=yval;
      }
   }

   // Set the fixed parameters again... this time they stick.
   transform1->SetFixedParameters(fparam);
   transform2->SetFixedParameters(fparam);
   std::cout << "Fixed parameters 1: " << transform1->GetFixedParameters() << std::endl;
   std::cout << "Fixed parameters 2: " << transform2->GetFixedParameters() << std::endl;

   transform1->SetParameters(param);
   transform2->SetParameters(param);
   std::cout << "Parameters 1: " << transform1->GetParameters() << std::endl;
   std::cout << "Parameters 2: " << transform2->GetParameters() << std::endl;

   std::string line;
   while(true) {
      std::getline(std::cin, line);
      std::stringstream ss(line);
      TransformType1::InputPointType p, q1,q2;
      ss >> p;
      if(!ss) break;
      std::cout << "Computing with " << transform1->GetNameOfClass() << std::endl;
      q1=transform1->TransformPoint(p);
      std::cout << "Computing with " << transform2->GetNameOfClass() << std::endl;
      q2=transform2->TransformPoint(p);
      std::cout << "Point is " << p << " transformed is " << q1 << " : " << q2 << " ( " << q2-q1 << " ) " << std::endl;
   }



   } catch (itk::ExceptionObject & e) {
      std::cerr << "Epic fail\n" << e << std::endl;
      return -2;
   }
   return EXIT_SUCCESS;
}

