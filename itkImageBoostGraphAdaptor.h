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

    OutputType Evaluate( const ImageType* Image,
                         const IndexType& iA,
                         const IndexType& iB ) const
      {
      PixelType a = Image->GetPixel( iA );
      PixelType b = Image->GetPixel( iB );

      return ( a -b ) * ( a -b );
      }

  };

template< class TInputImage,
          class TGraph,
          class TMetric > // Metric< TInputImage, typename TGraph::edge_property_type::value_type >
class ImageBoostGraphAdaptorBase : public LightObject
  {
public:
  typedef ImageBoostGraphAdaptorBase      Self;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;
  typedef LightObject                 Superclass;

  itkTypeMacro( ImageBoostGraphAdaptorBase, LightObject );

  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::SizeType       InputImageSizeType;
  typedef typename InputImageType::SizeValueType  InputImageSizeValueType;
  typedef typename InputImageType::IndexType      InputIndexType;

  typedef ConstShapedNeighborhoodIterator< InputImageType > NeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::OffsetType     NeighborhoodIteratorOffsetType;

  typedef TMetric MetricType;

  typedef TGraph GraphType;

  typedef typename GraphType::directed_selector     GraphDirectedType;

  typedef boost::graph_traits< GraphType >          GraphTraits;
  typedef typename GraphTraits::vertex_descriptor   VertexDescriptorType;
  typedef typename GraphTraits::edge_descriptor     EdgeDescriptorType;

  typedef typename GraphType::vertex_property_type  VertexPropertyType;
  typedef typename VertexPropertyType::value_type   VertexValueType;

  typedef typename GraphType::edge_property_type    EdgePropertyType;
  typedef typename EdgePropertyType::value_type     EdgeValueType;

  typedef typename boost::property_map< GraphType,
                                        boost::edge_weight_t >::type  WeightMapType;

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
    typename T::const_iterator end    = iOffsets.end();

    this->SetNeighbors( begin, end );
    }

  template< class TIterator >
  void SetNeighbors( const TIterator& iBegin, const TIterator& iEnd )
    {
    TIterator it = iBegin;

    while( it != iEnd )
      {
      this->m_OffsetList.push_back( *it );
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
  ImageBoostGraphAdaptorBase(){}
  virtual ~ImageBoostGraphAdaptorBase() {}

  GraphType               m_Graph;
  InputImageConstPointer  m_Image;
  MetricType              m_Metric;

  typedef std::list< NeighborhoodIteratorOffsetType > NeighborhoodIteratorOffsetContainerType;
  NeighborhoodIteratorOffsetContainerType m_OffsetList;

  NeighborhoodIteratorType GenerateNeighborhoodIterator( InputImageRegionType& oRegion )
    {
    if( !this->m_Image )
      {
      itkGenericExceptionMacro( << "input is null" );
      }

    oRegion = this->m_Image->GetRequestedRegion();

    InputImageSizeValueType numberOfVertices = oRegion.GetNumberOfPixels();
    this->m_Graph = GraphType( numberOfVertices );

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

    NeighborhoodIteratorType neighIt( radius, this->m_Image, oRegion );

    NeighborhoodIteratorOffsetType zeroOffset;
    zeroOffset.Fill( 0 );

    neighIt.DeactivateOffset( zeroOffset );

    for( typename NeighborhoodIteratorOffsetContainerType::const_iterator it = m_OffsetList.begin();
         it != m_OffsetList.end(); ++it )
      {
      neighIt.ActivateOffset( *it );
      }

    return neighIt;
    }

  virtual void GenerateData() = 0;

private:
  ImageBoostGraphAdaptorBase( const Self& );
  void operator = ( const Self& );
};

template< class TInputImage,
          class TOutEdgeListS,
          class TVertexListS,
          class TDirectedS,
          class TVertexProperty,
          class EdgeProperty,
          class GraphProperty,
          class EdgeListS,
          class TMetric >
class ImageBoostGraphAdaptor
  {};



template< class TInputImage,
          class TOutEdgeListS,
          class TVertexListS,
          class TVertexProperty,
          class TEdgeProperty,
          class TGraphProperty,
          class TEdgeListS,
          class TMetric > // Metric< TInputImage, typename TGraph::edge_property_type::value_type >
class ImageBoostGraphAdaptor<
    TInputImage,
    TOutEdgeListS,
    TVertexListS,
    boost::undirectedS,
    TVertexProperty,
    TEdgeProperty,
    TGraphProperty,
    TEdgeListS,
    TMetric >
    :
  public ImageBoostGraphAdaptorBase<
    TInputImage,
    boost::adjacency_list< TOutEdgeListS, TVertexListS, boost::undirectedS,
                           TVertexProperty, TEdgeProperty, TGraphProperty,
                           TEdgeListS >,
    TMetric
    >
{
public:
  typedef boost::adjacency_list< TOutEdgeListS, TVertexListS, boost::undirectedS,
    TVertexProperty, TEdgeProperty, TGraphProperty,
    TEdgeListS >                                                        GraphType;

  typedef ImageBoostGraphAdaptor                                        Self;
  typedef SmartPointer< Self >                                          Pointer;
  typedef SmartPointer< const Self >                                    ConstPointer;
  typedef ImageBoostGraphAdaptorBase< TInputImage, GraphType, TMetric > Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  itkTypeMacro( ImageBoostGraphAdaptor, ImageBoostGraphAdaptorBase );

  typedef typename Superclass::InputImageType           InputImageType;
  typedef typename Superclass::InputImageConstPointer   InputImageConstPointer;
  typedef typename Superclass::InputImageRegionType     InputImageRegionType;
  typedef typename Superclass::InputImageSizeType       InputImageSizeType;
  typedef typename Superclass::InputImageSizeValueType  InputImageSizeValueType;
  typedef typename Superclass::InputIndexType           InputIndexType;

  typedef typename Superclass::NeighborhoodIteratorType       NeighborhoodIteratorType;
  typedef typename Superclass::NeighborhoodIteratorOffsetType NeighborhoodIteratorOffsetType;

  typedef typename Superclass::MetricType MetricType;

  typedef typename Superclass::GraphDirectedType    GraphDirectedType;

  typedef typename Superclass::GraphTraits          GraphTraits;
  typedef typename Superclass::VertexDescriptorType VertexDescriptorType;
  typedef typename Superclass::EdgeDescriptorType   EdgeDescriptorType;

  typedef typename Superclass::VertexPropertyType VertexPropertyType;
  typedef typename Superclass::VertexValueType    VertexValueType;

  typedef typename Superclass::EdgePropertyType   EdgePropertyType;
  typedef typename Superclass::EdgeValueType      EdgeValueType;

  typedef typename Superclass::WeightMapType      WeightMapType;

protected:
  ImageBoostGraphAdaptor() {}
  ~ImageBoostGraphAdaptor() {}

  void GenerateData()
    {
    WeightMapType weightmap = get( boost::edge_weight, this->m_Graph );

    InputImageRegionType region;

    NeighborhoodIteratorType neighIt = this->GenerateNeighborhoodIterator( region );

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
//          itkAssertInDebugAndIgnoreInReleaseMacro( inside );

          EdgeDescriptorType e;

          std::pair< EdgeDescriptorType, bool> retrievedEdge = edge( u, v, this->m_Graph );

          if( !retrievedEdge.second )
            {
            bool inserted = false;
            boost::tie(e, inserted) = add_edge( u, v, this->m_Graph );
            weightmap[ e ] = this->m_Metric.Evaluate( this->m_Image, index, neighIndex );
            }
          }
        }
      }
    }

private:
  ImageBoostGraphAdaptor( const Self& );
  void operator = ( const Self& );
};


template< class TInputImage,
          class TOutEdgeListS,
          class TVertexListS,
          class TVertexProperty,
          class TEdgeProperty,
          class TGraphProperty,
          class TEdgeListS,
          class TMetric > // Metric< TInputImage, typename TGraph::edge_property_type::value_type >
class ImageBoostGraphAdaptor<
    TInputImage,
    TOutEdgeListS,
    TVertexListS,
    boost::directedS,
    TVertexProperty,
    TEdgeProperty,
    TGraphProperty,
    TEdgeListS,
    TMetric > :
  public ImageBoostGraphAdaptorBase<
    TInputImage,
    boost::adjacency_list< TOutEdgeListS, TVertexListS, boost::directedS,
                           TVertexProperty, TEdgeProperty, TGraphProperty,
                           TEdgeListS >,
    TMetric
    >
{
public:
  typedef boost::adjacency_list< TOutEdgeListS, TVertexListS, boost::directedS,
    TVertexProperty, TEdgeProperty, TGraphProperty,
    TEdgeListS >                                                        GraphType;

  typedef ImageBoostGraphAdaptor                                        Self;
  typedef SmartPointer< Self >                                          Pointer;
  typedef SmartPointer< const Self >                                    ConstPointer;
  typedef ImageBoostGraphAdaptorBase< TInputImage, GraphType, TMetric > Superclass;

  /** Method for creation through object factory */
  itkNewMacro( Self );

  itkTypeMacro( ImageBoostGraphAdaptor, ImageBoostGraphAdaptorBase );

  typedef typename Superclass::InputImageType           InputImageType;
  typedef typename Superclass::InputImageConstPointer   InputImageConstPointer;
  typedef typename Superclass::InputImageRegionType     InputImageRegionType;
  typedef typename Superclass::InputImageSizeType       InputImageSizeType;
  typedef typename Superclass::InputImageSizeValueType  InputImageSizeValueType;
  typedef typename Superclass::InputIndexType           InputIndexType;

  typedef typename Superclass::NeighborhoodIteratorType       NeighborhoodIteratorType;
  typedef typename Superclass::NeighborhoodIteratorOffsetType NeighborhoodIteratorOffsetType;

  typedef typename Superclass::MetricType MetricType;

  typedef typename Superclass::GraphDirectedType    GraphDirectedType;

  typedef typename Superclass::GraphTraits          GraphTraits;
  typedef typename Superclass::VertexDescriptorType VertexDescriptorType;
  typedef typename Superclass::EdgeDescriptorType   EdgeDescriptorType;

  typedef typename Superclass::VertexPropertyType VertexPropertyType;
  typedef typename Superclass::VertexValueType    VertexValueType;

  typedef typename Superclass::EdgePropertyType   EdgePropertyType;
  typedef typename Superclass::EdgeValueType      EdgeValueType;

  typedef typename Superclass::WeightMapType      WeightMapType;

protected:
  ImageBoostGraphAdaptor() {}
  ~ImageBoostGraphAdaptor() {}

  void GenerateData()
    {
    WeightMapType weightmap = get( boost::edge_weight, this->m_Graph );

    InputImageRegionType region;

    NeighborhoodIteratorType neighIt = this->GenerateNeighborhoodIterator( region );

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

          EdgeDescriptorType e;

          bool inserted = false;
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

