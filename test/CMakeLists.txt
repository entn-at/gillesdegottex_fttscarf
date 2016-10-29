cmake_minimum_required (VERSION 2.8)
project (tests)

# Add Boost library
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_definitions(-DBOOST_LOG_DYN_LINK)

# Add Voimoo library
include_directories(../include)
link_directories(${CMAKE_CURRENT_BINARY_DIR})

message("-- FFTScarf options")
# TODO Replace by fftscarf.h of the lib ?
if(FFTSCARF_PRECISION_FLOAT32)
    message("-- Single precision (32b)")
    add_definitions(-DFFTSCARF_PRECISION_FLOAT)
    add_definitions(-DOOFLOAT_FLOAT)
else(FFTSCARF_PRECISION_FLOAT32)
    message("-- Double precision (64b)")
endif(FFTSCARF_PRECISION_FLOAT32)

add_executable(test_version test_version.cpp)
target_link_libraries(test_version fftscarf ${Boost_LIBRARIES})

# add_executable(test_fft test_fft.cpp)
# #target_link_libraries(test_version ${Boost_LIBRARIES})
# target_link_libraries(test_fft fftscarf ${Boost_LIBRARIES})

# Optimization -----------------------------------------------------------------

# Speed up compilation
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

# -O2 -O3 -Ofast
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
