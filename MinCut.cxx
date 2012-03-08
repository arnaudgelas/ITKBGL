#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageBoostGraphAdaptor.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <boost/graph/one_bit_color_map.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>

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

  typedef double                                                              WeightType;

  typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property< boost::edge_weight_t, WeightType > > GraphType;

  typedef itk::IndexMetric< ImageType, WeightType >                           MetricType;
  typedef itk::ImageBoostGraphAdaptor< ImageType, GraphType, MetricType >     AdaptorType;

  std::vector< AdaptorType::NeighborhoodIteratorOffsetType > offset( 4 );

  size_t k = 0;
  offset[k][0] = -1;
  offset[k][1] = 0;
  k++;

  offset[k][0] = 0;
  offset[k][1] = -1;
  k++;

  offset[k][0] = 0;
  offset[k][1] = 1;
  k++;

  offset[k][0] = 1;
  offset[k][1] = 0;
  k++;

  AdaptorType::Pointer adaptor = AdaptorType::New();
  adaptor->SetInput( input );
  adaptor->SetNeighbors( offset );
  adaptor->Update();

  std::cout << "Graph constructed" << std::endl;

  typedef AdaptorType::GraphType GraphType;
  GraphType graph = adaptor->GetOutput();

  typedef AdaptorType::VertexDescriptorType   VertexDescriptorType;
  typedef AdaptorType::EdgeValueType          WeightType;

  // define a property map, `parities`, that will store a boolean value for each vertex.
  // Vertices that have the same parity after `stoer_wagner_min_cut` runs are on the same side of the min-cut.
  BOOST_AUTO( parities, boost::make_one_bit_color_map( num_vertices( graph ), get( boost::vertex_index, graph ) ) );

  // run the Stoer-Wagner algorithm to obtain the min-cut weight. `parities` is also filled in.
  int w = boost::stoer_wagner_min_cut( graph, get( boost::edge_weight, graph ), boost::parity_map( parities ) );

  std::cout << "The min-cut weight of G is " << w << ".\n" << std::endl;

  /*
  std::cout << "One set of vertices consists of:" << std::endl;
  size_t i;
  for( i = 0; i < num_vertices(g); ++i)
    {
    if( get(parities, i) )
      {
      std::cout << i << std::endl;
      }
    }

  std::cout << "The other set of vertices consists of:" << endl;
  for (i = 0; i < num_vertices(g); ++i)
    {
    if (!get(parities, i))
      {
      std::cout << i << std::endl;
      }
    }
  std::cout << std::endl;*/

  return EXIT_SUCCESS;
}
