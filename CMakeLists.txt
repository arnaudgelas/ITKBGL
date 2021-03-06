cmake_minimum_required( VERSION 2.8 )

project( ITKBGL )

find_package( ITK )
if( ${ITK_FOUND} )
  include( ${ITK_USE_FILE} )
endif()

find_package( Boost COMPONENTS graph )
include_directories( ${Boost_INCLUDE_DIRS} )

set(Boost_USE_MULTITHREADED ON)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ITKBGL_BINARY_DIR}/bin )

include( CTest )
add_executable( header header.cxx )
add_test( header ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/header )

add_executable( Undirected Undirected.cxx )
target_link_libraries( Undirected ${ITK_LIBRARIES} )

add_test( Undirected
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Undirected
  ${ITKBGL_SOURCE_DIR}/Data/Gourds.png
)

add_executable( Directed Directed.cxx )
target_link_libraries( Directed ${ITK_LIBRARIES} )

add_test( Directed
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Directed
  ${ITKBGL_SOURCE_DIR}/Data/Gourds.png
)

add_executable( Bidirectional Bidirectional.cxx )
target_link_libraries( Bidirectional ${ITK_LIBRARIES} )

add_test( Bidirectional
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Bidirectional
  ${ITKBGL_SOURCE_DIR}/Data/Gourds.png
)

add_executable( ShortestPath ShortestPath.cxx )
target_link_libraries( ShortestPath ${ITK_LIBRARIES} )

add_test( ShortestPath
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ShortestPath
  ${ITKBGL_SOURCE_DIR}/Data/Gourds.png
)

add_executable( MinCut MinCut.cxx )
target_link_libraries( MinCut ${ITK_LIBRARIES} )

#add_test( MinCut
#  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/MinCut
#  ${ITKBGL_SOURCE_DIR}/Data/Yinyang.png
#)

#add_executable( KolmogorovMaxFlow KolmogorovMaxFlow.cxx )
#target_link_libraries( KolmogorovMaxFlow ${ITK_LIBRARIES} )
#
#add_test( KolmogorovMaxFlow
#  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/KolmogorovMaxFlow
#  ${ITKBGL_SOURCE_DIR}/Data/Gourds.png
#)


