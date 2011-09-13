#include "itkGDCMSeriesFileNames.h"

int main( int ac, char* av[] )
{

  if ( ac != 2 )
    {
    std::cerr << "Usage: " << av[0] << " dicomDirectory ";
    return EXIT_FAILURE;
    }

  // Get the GDCM filenames from the directory
  itk::GDCMSeriesFileNames::Pointer names = itk::GDCMSeriesFileNames::New();
  names->RecursiveOff();
  names->SetUseSeriesDetails( true );
  names->SetLoadSequences( true );
  names->SetDirectory( av[1] );
  names->Update();

  std::vector<std::string> fileNames = names->GetInputFileNames();

  for( unsigned int d = 0; d < fileNames.size(); d++ )
    {
    std::cout << fileNames[d] << std::endl;
    }

  return EXIT_SUCCESS;

}
