#include <iostream>
#include <itkImage.h>
#include <itkBSplineDeformableTransformInitializer.h>
#include <itkBSplineDeformableTransform.h>

#define PATCHEDMETHOD

template<class ImageType>
typename ImageType::PointType GetImageCorner(typename ImageType::Pointer image, unsigned int c, int step=0) {
   const unsigned int dim=ImageType::ImageDimension;
   typename ImageType::SizeType sz=image->GetLargestPossibleRegion().GetSize();
   typename ImageType::IndexType cornerIndex;
   int pow2=1;
   for(unsigned int d=0;d<dim;++d) {
      cornerIndex[d]=((c & pow2)>0)?(sz[d]-1-step):0+step;
      pow2*=2;
   }
   typename ImageType::PointType p;
   image->TransformIndexToPhysicalPoint(cornerIndex, p);
   return p;
}


template<unsigned int dim, unsigned int order>
bool TestIt() {
   std::cout << "A BSpline test with dimension " << dim << " order " << order << std::endl;
   typedef itk::Image<unsigned char, dim> ImageType;
   typedef itk::BSplineDeformableTransform<double, dim, order> TransformType;
   typedef typename ImageType::SizeType SizeType;
   typedef typename ImageType::IndexType IndexType;
   typedef typename ImageType::PointType PointType;
   typedef typename ImageType::SpacingType SpacingType;
   typedef typename ImageType::RegionType RegionType;
   typedef typename ImageType::DirectionType DirectionType;
   typedef typename TransformType::ParametersType ParametersType;

   typedef itk::BSplineDeformableTransformInitializer<TransformType,ImageType> TransformInitializerType;

   typename ImageType::Pointer anImage=ImageType::New();
   {
      SizeType   aSize;
      PointType  anOrigin;
      SpacingType aSpacing;
      for(unsigned int i=0;i<dim;++i) {aSize[i]=50+25*i; anOrigin[i]=-11.1+5*i+double(i)/8;aSpacing=0.75+double(i)/8;}
      // avoid zeros or equal sizes and spacings, just to be sure no one mixed up...

      DirectionType aDirection;
      aDirection.SetIdentity();
      //Rotate a little around each dir

      IndexType aStart;aStart.Fill(0);
      RegionType aRegion(aStart,aSize);
      anImage->SetRegions(aRegion);
      anImage->SetOrigin(anOrigin);
      anImage->SetSpacing(aSpacing);
      anImage->SetDirection(aDirection);
      anImage->Allocate();
   }

   std::cout << "The image given to the BSpline is: \n" << anImage << std::endl;

   typename TransformType::Pointer aTransform=TransformType::New();

   typename TransformInitializerType::Pointer anInitializer=TransformInitializerType::New();
   anInitializer->SetImage(anImage);
   anInitializer->SetTransform(aTransform);
   typename TransformInitializerType::TransformType::SizeType aTransformSize;
   for(unsigned int i=0;i<dim;++i) {aTransformSize[i]=3;}//+i;}
#ifdef PATCHEDMETHOD
   aTransform->SetTransformDomainMeshSize(aTransformSize);
#else
   anInitializer->SetGridSizeInsideTheImage(aTransformSize);
#endif
   anInitializer->InitializeTransform();

   const unsigned int NP=aTransform->GetNumberOfParameters();
   ParametersType p=aTransform->GetParameters();
   assert(p.Size()==NP);
   std::cout << "Number of Transform Parameters: " << NP << std::endl;

   const unsigned int NPperDim=NP/dim;
   assert(NPperDim*dim==NP);

   for(unsigned int np=0;np<NPperDim;++np) {
      for(unsigned int d=0;d<dim;++d) {
         const unsigned int index=np+d*NPperDim;
         assert (index<NP);
         p[index]=5+double(np)/NPperDim;
      }
   }

   aTransform->SetParameters(p);

   ParametersType f=aTransform->GetFixedParameters();
   std::cout << "Number of Transform Parameters: " << f.Size() << std::endl;
   std::cout << "Fixed Parameters: " << f << std::endl;

#if ITK_VERSION_MAJOR == 3
   std::cout << "Debug: ITK3\n";
 //this doesnt work in ITK3.20
   typename TransformType::ImagePointer * coeff=aTransform->GetCoefficientImage();
#else
   std::cout << "Debug: ITK4\n";
#ifdef PATCHEDMETHOD
   typename TransformType::CoefficientImageArray coeff=aTransform->GetCoefficientImages();
#else
   typename TransformType::CoefficientImageArray coeff=aTransform->GetCoefficientImage();
#endif
#endif
   // Dump the positions of the grid points
   for(unsigned int i=0;i<dim;++i) {
      std::cout << "Positions of grid points along dimension " << i << std::endl;
      PointType p1=GetImageCorner<ImageType>(anImage,static_cast<unsigned int>(0),false);
      PointType p2=GetImageCorner<ImageType>(anImage,
        static_cast<unsigned int>(std::pow(2.,static_cast<int>( i ))), false );
      std::cout << "ImageLimits: " << p1[i] << " -> " << p2[i] << std::endl;
      typename TransformType::ImageType::Pointer im = coeff[i];
      assert(im.IsNotNull());
      typename TransformType::ImageType::SizeType sz=im->GetLargestPossibleRegion().GetSize();
      typename TransformType::ImageType::IndexType ix; ix.Fill(0);
      for(unsigned int j=0; j<sz[i];++j) {
         ix[i]=j;
         PointType p;
         im->TransformIndexToPhysicalPoint(ix,p);
         std::cout << "    " << j << ": " << p;
      }
      std::cout << std::endl;
   }


   // Check the image corners - are they in or out.
   unsigned int c=0;
   unsigned int nbCorners=std::pow(2.,static_cast<int>( dim ));
   while(c<nbCorners) {
      SizeType sz=anImage->GetLargestPossibleRegion().GetSize();
//      PointType pin=GetImageCorner<ImageType>(anImage,c,1);
//      std::cout << "Origin: " << anImage->GetOrigin() << " point " << pin << std::endl;
//      PointType tin=aTransform->TransformPoint(pin);
//      {
//         std::cout << "Corner " << c << " pin " << pin << " " << tin << std::endl;
//      }
      PointType pex=GetImageCorner<ImageType>(anImage,c,0);
      PointType tex=aTransform->TransformPoint(pex);
      {
         std::cout << "Corner " << c << " pex " << pex << " " << tex << std::endl;
      }
//      PointType pout=GetImageCorner<ImageType>(anImage,c,-1);
//      PointType tout=aTransform->TransformPoint(pout);
//      {
//         std::cout << "Corner " << c << " pout " << pout << " " << tout << std::endl;
//      }
      c=c+1;
   }

   return true;

}


int main(int argc, char ** argv) {
//   TestIt<2,1>();
//   TestIt<2,2>();
   TestIt<2,3>();
   //TestIt<2,4>();
   return EXIT_SUCCESS;
}
