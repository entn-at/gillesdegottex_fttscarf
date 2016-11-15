cmake_minimum_required (VERSION 3.2)
project (test)

# Add Boost unit test library.
set(Boost_USE_STATIC_LIBS on)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Add FFTScarf library
include_directories(../include)
link_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_version test_version.cpp)
target_link_libraries(test_version fftscarf ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_executable(test_fftlibs test_fftlibs.cpp)
target_link_libraries(test_fftlibs fftscarf ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})


# Add support for unit tests e.g. ctest.
enable_testing()
add_test(test_version test_version)
add_test(test_fftlibs test_fftlibs)

# Optimization -----------------------------------------------------------------

# Speed up compilation
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

# -O2 -O3 -Ofast
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
