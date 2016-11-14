cmake_minimum_required (VERSION 2.8)
project (tests)

# Add Boost library
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_definitions(-DBOOST_LOG_DYN_LINK)
#set(Boost_USE_STATIC_LIBS on)

# Add FFTScarf library
include_directories(../include)
link_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_version test_version.cpp)
target_link_libraries(test_version fftscarf ${Boost_LIBRARIES})

add_executable(test_fftlibs test_fftlibs.cpp)
target_link_libraries(test_fftlibs fftscarf ${Boost_LIBRARIES})

# Optimization -----------------------------------------------------------------

# Speed up compilation
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

# -O2 -O3 -Ofast
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
