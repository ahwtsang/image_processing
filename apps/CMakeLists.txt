
project (apps)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)

find_package(ITK 5.2.1 EXACT REQUIRED)
include(${ITK_USE_FILE})

add_subdirectory(ImageOrienter)
add_subdirectory(N4BiasCorrection)