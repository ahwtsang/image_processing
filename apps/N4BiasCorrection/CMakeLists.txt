
project(N4BiasCorrection)

SET(N4BiasFieldCorrectionApp_SRCS
  N4BiasCorrection.cpp
)

ADD_EXECUTABLE(N4BiasFieldCorrectionApp ${N4BiasFieldCorrectionApp_SRCS})

TARGET_LINK_LIBRARIES(N4BiasFieldCorrectionApp ${ITK_LIBRARIES})
