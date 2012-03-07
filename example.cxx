#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageBoostGraphAdaptor.h"
#include "itkImageRegionConstIteratorWithIndex.h"

int main( int argc, char* argv[] )
  {
  if( argc != 2 )
    {
    std::cerr << argv[0] << " <InputImage>" << std::endl;
    return EXIT_FAILURE;
    }
  typedef unsigned char PixelType;
  const unsigned int Dimension = 2;

  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::Pointer input = reader->GetOutput();

  typedef itk::ImageBoostGraphAdaptor< ImageType > AdaptorType;
  AdaptorType::Pointer adaptor = AdaptorType::New();
  adaptor->SetInput( input );
  adaptor->Update();

  typedef AdaptorType::GraphType GraphType;
  GraphType graph = adaptor->GetOutput();

  ImageType::RegionType region = input->GetLargestPossibleRegion();

  typedef AdaptorType::VertexDescriptorType VertexDescriptorType;

  typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;
  IteratorType it( input, region );
  it.GoToBegin();

  while( !it.IsAtEnd() )
    {
    ImageType::IndexType idx = it.GetIndex();
    bool inside = false;

    VertexDescriptorType v = adaptor->GetVertexFromIndex( idx, inside );

    if( inside )
      {
      ImageType::IndexType test = adaptor->GetIndexFromVertex( v );

      if( test != idx )
        {
        std::cerr << "test: " << test << std::endl;
        std::cerr << "idx: " << idx << std::endl;
        std::cerr << "v: " << v << std::endl;
        return EXIT_FAILURE;
        }
      }
    else
      {
      std::cerr << "Outside: " << idx << std::endl;
      return EXIT_FAILURE;
      }
    ++it;
    }

  std::cout << "SUCCESS!" << std::endl;
  return EXIT_SUCCESS;
  }
