#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageBoostGraphAdaptor.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

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

  std::vector< AdaptorType::NeighborhoodIteratorOffsetType > offset( 8 );

  size_t k = 0;
  offset[k][0] = -1;
  offset[k][1] = -1;
  k++;

  offset[k][0] = -1;
  offset[k][1] = 0;
  k++;

  offset[k][0] = -1;
  offset[k][1] = 1;
  k++;

  offset[k][0] = 0;
  offset[k][1] = -1;
  k++;

  offset[k][0] = 0;
  offset[k][1] = 1;
  k++;

  offset[k][0] = 1;
  offset[k][1] = -1;
  k++;

  offset[k][0] = 1;
  offset[k][1] = 0;
  k++;

  offset[k][0] = 1;
  offset[k][1] = 1;
  k++;


  AdaptorType::Pointer adaptor = AdaptorType::New();
  adaptor->SetInput( input );
  adaptor->SetNeighbors( offset );
  adaptor->Update();

  std::cout << "Graph constructed" << std::endl;

  typedef AdaptorType::GraphType GraphType;
  GraphType graph = adaptor->GetOutput();

  typedef AdaptorType::VertexDescriptorType   VertexDescriptorType;
  typedef AdaptorType::EdgeDescriptorType     EdgeDescriptorType;
  typedef AdaptorType::WeightType             WeightType;

  VertexDescriptorType s, t;

  flow = boykov_kolmogorov_max_flow(graph, s, t); // a list of sources will be returned in s, and a list of sinks will be returned in t

  return EXIT_SUCCESS;
}

