# Sources of non-templated classes.
project(ParametersParserLib)

set( param_SRCS
  itkParameterFileParser.h
  itkParameterFileParser.cxx
  itkParameterMapInterface.h
  itkParameterMapInterface.cxx
)

## source_group( "Parser" FILES ${param_SRCS} )

add_library( ParametersParser STATIC ${param_SRCS} )

target_link_libraries( ParametersParser ${ITK_LIBRARIES} )
