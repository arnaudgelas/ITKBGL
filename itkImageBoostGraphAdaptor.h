#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "itkLightObject.h"
#include "itkConstShapedNeighborhoodIterator.h"

namespace itk
{
template< class TImage, class TOutput >
class IndexMetric
  {
public:
    typedef TImage ImageType;
    typedef TOutput OutputType;

    typedef typename ImageType::IndexType IndexType;
    typedef typename ImageType::PixelType PixelType;

    OutputType Evaluate( const ImageType* Image, const IndexType& iA, const IndexType& iB ) const
      {
      PixelType a = Image->GetPixel( iA );
      PixelType b = Image->GetPixel( iB );

      return ( a -b ) * ( a -b );
      }

  };

template< class TInputImage >
class ImageBoostGraphAdaptor : public LightObject
  {
public:
  typedef ImageBoostGraphAdaptor      Self;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;
  typedef LightObject                 Superclass;

  itkNewMacro( Self );
  itkTypeMacro( ImageBoostGraphAdaptor, LightObject );

  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::SizeType       InputImageSizeType;
  typedef typename InputImageType::SizeValueType  InputImageSizeValueType;
  typedef typename InputImageType::IndexType      InputIndexType;

  typedef ConstShapedNeighborhoodIterator< InputImageType > NeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::OffsetType     NeighborhoodIteratorOffsetType;

  typedef double WeightType;

  typedef IndexMetric< TInputImage, WeightType > MetricType;

  typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property, boost::property< boost::edge_weight_t, WeightType > > GraphType;

  typedef boost::graph_traits< GraphType >                                GraphTraits;
  typedef typename GraphTraits::vertex_descriptor                         VertexDescriptorType;
  typedef typename GraphTraits::edge_descriptor                           EdgeDescriptorType;
  typedef typename boost::property_map< GraphType,
                                        boost::edge_weight_t >::type      WeightMapType;

  void SetInput( const InputImageType* Image )
    {
    this->m_Image = Image;
    }

  const InputImageType* GetInput()
    {
    return static_cast< const InputImageType* >( this->m_Image.GetPointer() );
    }

  template< class T >
  void SetNeighbors( const T & iOffsets )
    {
    typename T::const_iterator begin  = iOffsets.begin();
    typename T::const_iterator end    = iOffsets.begin();

    this->SetNeighbors( begin, end );
    }

  template< class TIterator >
  void SetNeighbors( const TIterator& iBegin, const TIterator& iEnd )
    {
    TIterator it = iBegin;

    while( it != iEnd )
      {
      this->m_OffsetList.insert( *it );
      ++it;
      }
    }

  void Update()
    {
    this->GenerateData();
    }

  const GraphType & GetOutput() const
    {
    return this->m_Graph;
    }

  VertexDescriptorType GetVertexFromIndex( const InputIndexType& idx,
                                           bool& oIsInside ) const
    {
    typename InputImageType::OffsetValueType res = 0;
    InputImageRegionType region = this->m_Image->GetLargestPossibleRegion();
    oIsInside = region.IsInside( idx );
    if( oIsInside )
      {
      res = this->m_Image->ComputeOffset( idx );
      }
    return vertex( res, this->m_Graph );
    }

  InputIndexType GetIndexFromVertex( const VertexDescriptorType& iV )
    {
    return this->m_Image->ComputeIndex( static_cast< typename InputImageType::OffsetValueType >( iV ) );
    }

protected:
  ImageBoostGraphAdaptor(){}
  ~ImageBoostGraphAdaptor() {}

  GraphType               m_Graph;
  InputImageConstPointer  m_Image;
  InputImageSizeType      m_Size;
  MetricType              m_Metric;

  typedef std::set< NeighborhoodIteratorOffsetType > NeighborhoodIteratorOffsetContainerType;
  NeighborhoodIteratorOffsetContainerType m_OffsetList;

  void GenerateData()
  {
    if( !this->m_Image )
      {
      itkGenericExceptionMacro( << "input is null" );
      }

    InputImageRegionType region = this->m_Image->GetRequestedRegion();

    m_Size = region.GetSize();

    InputImageSizeValueType numberOfVertices = region.GetNumberOfPixels();
    this->m_Graph = GraphType( numberOfVertices );

    WeightMapType weightmap = get( boost::edge_weight, this->m_Graph );

    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 0 );

    for( typename NeighborhoodIteratorOffsetContainerType::const_iterator it = m_OffsetList.begin();
         it != m_OffsetList.end(); ++it )
      {
      for( unsigned int dim = 0; dim < InputImageType::ImageDimension; ++dim )
        {
        typename NeighborhoodIteratorOffsetType::OffsetValueType absValue = vnl_math_abs( ( *it )[ dim ] );
        if( absValue > radius[ dim ] )
          {
          radius[ dim ] = absValue;
          }
        }
      }

    NeighborhoodIteratorType neighIt( radius, this->m_Image, region );

    for( typename NeighborhoodIteratorOffsetContainerType::const_iterator it = m_OffsetList.begin();
         it != m_OffsetList.end(); ++it )
      {
      neighIt.ActivateOffset( *it );
      }

    for( neighIt.GoToBegin(); !neighIt.IsAtEnd(); ++neighIt )
      {
      InputIndexType index = neighIt.GetIndex();
      bool inside = false;

      VertexDescriptorType u = this->GetVertexFromIndex( index, inside );

//      itkAssertInDebugAndIgnoreInReleaseMacro( inside );

      for( typename NeighborhoodIteratorType::ConstIterator internalIt = neighIt.Begin();
           !internalIt.IsAtEnd(); ++internalIt )
        {
        unsigned int i = internalIt.GetNeighborhoodIndex();
        InputIndexType neighIndex = neighIt.GetIndex( i );

        bool IsInBounds = region.IsInside( neighIndex );

        if( IsInBounds && index != neighIndex )
          {
          VertexDescriptorType v = this->GetVertexFromIndex( neighIndex, inside );
          itkAssertInDebugAndIgnoreInReleaseMacro( inside );

          EdgeDescriptorType e;

          bool inserted;
          boost::tie(e, inserted) = add_edge( u, v, this->m_Graph );
          weightmap[ e ] = this->m_Metric.Evaluate( this->m_Image, index, neighIndex );
          }
        }
      }
  }

private:
  ImageBoostGraphAdaptor( const Self& );
  void operator = ( const Self& );
};

}

