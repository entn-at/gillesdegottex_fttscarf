cmake_minimum_required (VERSION 3.2)
project (test)

# Add Boost unit test library.
set(Boost_USE_STATIC_LIBS on)
find_package(Boost COMPONENTS unit_test_framework program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Add FFTScarf library
link_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_version test_version.cpp)
target_link_libraries(test_version fftscarf ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(test_version ${Boost_LIBRARIES})
if(WIN32)
    if(FFTSCARF_PRECISION_SINGLE)
        target_link_libraries(test_version C:/projects/fftscarf/lib/libfft/libfftw3f-3.lib)
    endif()
    if(FFTSCARF_PRECISION_DOUBLE)
        target_link_libraries(test_version C:/projects/fftscarf/lib/libfft/libfftw3-3.lib)
    endif()
endif()

add_executable(test_misc test_misc.cpp)
target_link_libraries(test_misc fftscarf ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(test_misc ${Boost_LIBRARIES})
if(WIN32)
    if(FFTSCARF_PRECISION_SINGLE)
        target_link_libraries(test_misc C:/projects/fftscarf/lib/libfft/libfftw3f-3.lib)
    endif()
    if(FFTSCARF_PRECISION_DOUBLE)
        target_link_libraries(test_misc C:/projects/fftscarf/lib/libfft/libfftw3-3.lib)
    endif()
endif()

add_executable(test_fftlibs test_fftlibs.cpp)
target_link_libraries(test_fftlibs fftscarf ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(test_fftlibs ${Boost_LIBRARIES})
if(WIN32)
    if(FFTSCARF_PRECISION_SINGLE)
        target_link_libraries(test_fftlibs C:/projects/fftscarf/lib/libfft/libfftw3f-3.lib)
    endif()
    if(FFTSCARF_PRECISION_DOUBLE)
        target_link_libraries(test_fftlibs C:/projects/fftscarf/lib/libfft/libfftw3-3.lib)
    endif()
endif()

# Add tests using ctest
enable_testing()
add_test(test_version test_version)
add_test(test_fftlibs test_misc)
add_test(test_fftlibs test_fftlibs)
#add_test(test_fftlibs32 test_fftlibs --size=32) // Doesn't run on OSX ...
#add_test(test_fftlibs4096 test_fftlibs --size=4096) // Doesn't run on OSX ...

# Optimization -----------------------------------------------------------------

if(WIN32) # TODO Should be test on GCC because the followin might work on mingw
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /O2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2")
else() # TODO Should be test on GCC because the followin might work on mingw
    # Speed up compilation
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

    # -O2 -O3 -Ofast
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()
