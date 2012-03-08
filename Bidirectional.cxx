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

  typedef double                                                              WeightType;

  typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS,
    boost::no_property, boost::property< boost::edge_weight_t, WeightType > > GraphType;

  typedef itk::IndexMetric< ImageType, WeightType >                           MetricType;
  typedef itk::ImageBoostGraphAdaptor< ImageType, GraphType, MetricType >     AdaptorType;
  AdaptorType::Pointer adaptor = AdaptorType::New();
  adaptor->SetInput( input );

  std::vector< AdaptorType::NeighborhoodIteratorOffsetType > offset( 2 );

  size_t k = 0;
  offset[k][0] = -1;
  offset[k][1] = 0;
  k++;

  offset[k][0] = 1;
  offset[k][1] = 0;
  k++;

  adaptor->SetNeighbors( offset );
  adaptor->Update();

  typedef AdaptorType::GraphType GraphType;
  GraphType graph = adaptor->GetOutput();

  ImageType::RegionType region = input->GetLargestPossibleRegion();

  ImageType::SizeType size = region.GetSize();

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

    int oDeg = boost::out_degree( v, graph );
    if( oDeg != 2 )
      {
      if( oDeg == 1 )
        {
        if( ( idx[0] != 0 ) && ( idx[1] == 0 ) && ( idx[0] != size[0] -1 ) && ( idx[1] != size[1] - 1 ) )
          {
          std::cout << "degree[ " << idx << " ] = " << oDeg << std::endl;
          return EXIT_FAILURE;
          }
        }
      else
        {
        std::cout << "degree[ " << idx << " ] = " << oDeg << std::endl;
        return EXIT_FAILURE;
        }
      }

    int iDeg = boost::in_degree( v, graph );
    if( iDeg != 2 )
      {
      if( iDeg == 1 )
        {
        if( ( idx[0] != 0 ) && ( idx[1] == 0 ) && ( idx[0] != size[0] -1 ) && ( idx[1] != size[1] - 1 ) )
          {
          std::cout << "degree[ " << idx << " ] = " << iDeg << std::endl;
          return EXIT_FAILURE;
          }
        }
      else
        {
        std::cout << "degree[ " << idx << " ] = " << iDeg << std::endl;
        return EXIT_FAILURE;
        }
      }

    bool found = false;
    AdaptorType::EdgeDescriptorType e;
    boost::tie( e, found ) = boost::edge( v, v, graph );

    if( found )
      {
      std::cout << v << std::endl;
      return EXIT_FAILURE;
      }
    ++it;
    }



  std::cout << "SUCCESS!" << std::endl;
  return EXIT_SUCCESS;
  }

