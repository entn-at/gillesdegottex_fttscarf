cmake_minimum_required (VERSION 2.8)
project (fftscarf)

add_library(fftscarf src/fftscarf.cpp)

message("\n-- FFTScarf configuring lib")

option(FFTSCARF_PRECISION_FLOAT32 "Precision" ON)
if(FFTSCARF_PRECISION_FLOAT32)
    message("--   Single precision (32b)")
#     add_definitions(-DFFTSCARF_PRECISION_FLOAT32)
else(FFTSCARF_PRECISION_FLOAT32)
    message("--   Double precision (64b)")
endif(FFTSCARF_PRECISION_FLOAT32)

message("--   FFT Implementation")
option(FFTSCARF_FFT_OOURA "Ooura's FFT" ON)
option(FFTSCARF_FFT_FFTREAL "FFTReal" ON)
option(FFTSCARF_FFT_PFFFT "PFFFT" ON)
option(FFTSCARF_FFT_DJBFFT "DJBFFT" OFF)
option(FFTSCARF_FFT_FFTW3 "FFTW3" OFF)
option(FFTSCARF_FFT_FFTS "FFTS" ON)
option(FFTSCARF_FFT_IPP "IPP" OFF)

if(FFTSCARF_FFT_OOURA)
    message("--     Using Ooura (built-in)")
    if(FFTSCARF_PRECISION_FLOAT32)
        add_definitions(-DOOFLOAT_FLOAT)
    endif(FFTSCARF_PRECISION_FLOAT32)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_ooura.cpp)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ooura/fftsg.c)
#    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ooura/fft4g.cpp)
endif(FFTSCARF_FFT_OOURA)

if(FFTSCARF_FFT_FFTREAL)
    message("--     Using FFTReal (built-in)")
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_fftreal.cpp)
endif(FFTSCARF_FFT_FFTREAL)

if(FFTSCARF_FFT_PFFFT)
    message("--     Using Pretty Fast FFT (PFFFT) (built-in)")
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_pffft.cpp)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/pffft/pffft.c)
endif(FFTSCARF_FFT_PFFFT)

# if(FFTSCARF_FFT_DJBFFT)
#     message("--     Using DJBFFT (shipped with)")
# #    message("--     Configuring DJBFFT library")
#     target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_djbfft.cpp)
#     target_link_libraries(fftscarf ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/djbfft-0.76/djbfft.a)
# endif(FFTSCARF_FFT_DJBFFT)

if(FFTSCARF_FFT_FFTS)
    message("--     Using FFTS (shipped with)")
    message("--         Configuring FFTS library...")
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_ffts.cpp)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ffts ${CMAKE_CURRENT_BINARY_DIR}/ffts_build)
    add_custom_target(libffts
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ffts_build/
                     COMMAND make)
    add_dependencies(fftscarf libffts)
    target_link_libraries(fftscarf ${CMAKE_CURRENT_BINARY_DIR}/ffts_build/libffts.a)
endif(FFTSCARF_FFT_FFTS)

if(FFTSCARF_FFT_FFTW3)
    message("--     Using FFTW3 (external lib)")
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_fftw3.cpp)
    if(FFTSCARF_PRECISION_FLOAT32)
        target_link_libraries(fftscarf fftw3f)
    else(FFTSCARF_PRECISION_FLOAT32)
        target_link_libraries(fftscarf fftw3)
    endif(FFTSCARF_PRECISION_FLOAT32)
endif(FFTSCARF_FFT_FFTW3)

# if(FFTSCARF_FFT_IPP)
#     message("--     Using Intel Integrated Performance Primitives (Intel IPP) (external lib)")
#     include_directories(/opt/intel/ipp/include)
#     target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_ipp.cpp)
#     target_link_libraries(fftscarf /opt/intel/ipp/lib/intel64/libipps.a /opt/intel/ipp/lib/intel64/libippcore.a  /opt/intel/ipp/lib/intel64/libippvm.a)
# endif(FFTSCARF_FFT_IPP)

configure_file(include/fftscarf.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/fftscarf.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(include)


# Testing ---------------------------------------------------------------------

option(FFTSCARF_BUILD_TEST "Build tests" OFF)
if(FFTSCARF_BUILD_TEST)
    message("\n-- FFTScarf configuring tests")
    add_subdirectory(test)

    # Add a custom target for running all the tests
    add_custom_target(test_run)
    add_dependencies(test_run test_run_version test_run_fftlibs)

    add_custom_target(test_run_version
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test
                        COMMAND ./test_version)
    add_dependencies(test_run_version test_version)

     add_custom_target(test_run_fftlibs
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test
                         COMMAND ./test_fftlibs)
     add_dependencies(test_run_fftlibs test_fftlibs)

endif(FFTSCARF_BUILD_TEST)

message("")

# Optimization -----------------------------------------------------------------

# Speed up compilation
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

# -O2 -O3 -Ofast
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Install ---------------------------------------------------------------------
message("-- Installation path: " ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION include/
        FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include/
        FILES_MATCHING PATTERN "*.h")
install(TARGETS fftscarf DESTINATION lib)
