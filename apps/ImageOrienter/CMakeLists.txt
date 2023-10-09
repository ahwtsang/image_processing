project(ImageOrienter)

SET(ImageOrienterApp_SRCS
  ImageOrienter.cpp
)

ADD_EXECUTABLE(ImageOrienterApp ${ImageOrienterApp_SRCS})

###TARGET_INCLUDE_DIRECTORIES(ImageReorienterApp PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/libs/ParameterFileParserLib>)
include_directories(${CMAKE_SOURCE_DIR}/ParameterFileParserLib)

TARGET_LINK_LIBRARIES(ImageOrienterApp
 ParametersParser
 ${ITK_LIBRARIES}
)
