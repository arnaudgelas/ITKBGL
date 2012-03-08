#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageBoostGraphAdaptor.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <boost/graph/dijkstra_shortest_paths.hpp>

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

  ImageType::IndexType idx1, idx2;
  idx1[0] = 320;
  idx1[1] = 240;

  idx2[0] = 160;
  idx2[1] = 120;

  typedef AdaptorType::VertexDescriptorType   VertexDescriptorType;
  typedef AdaptorType::WeightType             WeightType;

  bool inside = false;
  VertexDescriptorType v1 = adaptor->GetVertexFromIndex( idx1, inside );
  VertexDescriptorType v2 = adaptor->GetVertexFromIndex( idx2, inside );

  typedef boost::property_map< GraphType, boost::vertex_index_t > IndexMapType;

  typedef boost::iterator_property_map< VertexDescriptorType*, IndexMapType, VertexDescriptorType, VertexDescriptorType& >  PredecessorMapType;
  typedef boost::iterator_property_map < WeightType*, IndexMapType, WeightType, WeightType& >                               DistanceMapType;

  std::vector< VertexDescriptorType > Predecessors( num_vertices( graph ) );
  std::vector< WeightType >           Distances( num_vertices( graph ) );

  boost::dijkstra_shortest_paths(graph, v1,
                                 boost::predecessor_map( &Predecessors[0] ).distance_map( &Distances[0] ) );


  typedef std::vector<boost::graph_traits< GraphType >::edge_descriptor> PathType;
  PathType path;

  VertexDescriptorType v = v2;                      // We want to start at the destination and work our way back to the source
  for( VertexDescriptorType u = Predecessors[ v ];  // Start by setting 'u' to the destintaion node's predecessor
       u != v;                                     // Keep tracking the path until we get to the source
       v = u, u = Predecessors[ v ] )               // Set the current vertex to the current predecessor, and the predecessor to one level up
    {
    std::pair<boost::graph_traits< GraphType >::edge_descriptor, bool> edgePair = edge( u, v, graph );
    path.push_back( edgePair.first );
    }

  // Write shortest path
  std::cout << "Shortest path:" << std::endl;
  double totalDistance = 0;
  for( PathType::reverse_iterator pathIterator = path.rbegin();
       pathIterator != path.rend();
       ++pathIterator )
    {
    VertexDescriptorType edge_source      = source( *pathIterator, graph );
    VertexDescriptorType edge_destination = source( *pathIterator, graph );

    std::cout << adaptor->GetIndexFromVertex( edge_source ) << " -> "
              << adaptor->GetIndexFromVertex( edge_destination )
              << " = " << get( boost::edge_weight, graph, *pathIterator ) << std::endl;
    }

  std::cout << std::endl;
  std::cout << "Distance: " << Distances[v2] << std::endl;


  return EXIT_SUCCESS;
}
