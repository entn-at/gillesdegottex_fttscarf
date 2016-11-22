cmake_minimum_required (VERSION 3.2)
project (fftscarf)

add_library(fftscarf src/fftscarf.cpp)

message("\n-- FFTScarf configuring lib")

option(FFTSCARF_PRECISION_SINGLE "Enable single precision" ON)
option(FFTSCARF_PRECISION_DOUBLE "Enable double precision" ON)
option(FFTSCARF_PRECISION_LONGDOUBLE "Enable long double precision" OFF)
option(FFTSCARF_PRECISION_DEFAULT "Set default precision (32, 64, 128)" 32)

option(FFTSCARF_PLAN_PROTECTACCESS "Protect access of the plan resize and execution using a mutex (need boost library)." OFF)
if(FFTSCARF_PLAN_PROTECTACCESS)
    message("-- Protect plan access (resize and execution) using a mutex (using Boost library)")
    set(Boost_USE_STATIC_LIBS on)
    find_package(Boost COMPONENTS thread REQUIRED)
    add_definitions(-DBOOST_DATE_TIME_POSIX_TIME_STD_CONFIG)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries(fftscarf ${Boost_LIBRARIES})
else()
    message("-- Plan access are not protected")
endif()

if(FFTSCARF_PRECISION_SINGLE)
    message("-- Enable single precision (float 32b)")
else()
    set(FFTSCARF_FFT_FFTS OFF CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_FFT_PFFFT OFF CACHE BOOL "docstring" FORCE)
endif()
if(FFTSCARF_PRECISION_DOUBLE)
    message("-- Enable double precision (float 64b)")
endif(FFTSCARF_PRECISION_DOUBLE)
if(FFTSCARF_PRECISION_LONGDOUBLE)
    message("-- Enable long double precision (float 128b)")
endif(FFTSCARF_PRECISION_LONGDOUBLE)
if(FFTSCARF_PRECISION_DEFAULT EQUAL 32)
    message("-- Default is single precision (float 32b)")
endif()
if(FFTSCARF_PRECISION_DEFAULT EQUAL 64)
    message("-- Default is double precision (float 64b)")
endif()
if(FFTSCARF_PRECISION_DEFAULT EQUAL 128)
    message("-- Default is long double precision (float 128b)")
endif()

option(FFTSCARF_FFT_IPP "IPP" OFF)
option(FFTSCARF_FFT_FFTS "FFTS" OFF)
option(FFTSCARF_FFT_PFFFT "PFFFT" OFF)
option(FFTSCARF_FFT_FFTW3 "FFTW3" OFF)
option(FFTSCARF_FFT_OOURA "Ooura's FFT" ON)
option(FFTSCARF_FFT_FFTREAL "FFTReal" OFF)
option(FFTSCARF_FFT_DFT "DFT" OFF)
#option(FFTSCARF_FFT_DJBFFT "DJBFFT" OFF)

if(FFTSCARF_BUILD_TEST)
    set(FFTSCARF_FFT_FFTW3 ON CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_PRECISION_LONGDOUBLE ON CACHE BOOL "docstring" FORCE)
endif()

# FFT Implementations --------------------------------------------------------

message("-- FFT Implementation")

function(cat IN_FILE OUT_FILE)
  file(READ ${IN_FILE} CONTENTS)
  file(APPEND ${OUT_FILE} "${CONTENTS}")
endfunction()

configure_file(include/fftscarf.h.in ${CMAKE_CURRENT_BINARY_DIR}/fftscarf.h)

if(FFTSCARF_FFT_IPP)
    message("--   Enable Intel Integrated Performance Primitives (Intel IPP) (external lib)")
    cat(${CMAKE_CURRENT_SOURCE_DIR}/include/fftscarf_ipp.h ${CMAKE_CURRENT_BINARY_DIR}/fftscarf.h)
    include_directories(/opt/intel/ipp/include)
    target_link_libraries(fftscarf /opt/intel/ipp/lib/intel64/libipps.a /opt/intel/ipp/lib/intel64/libippcore.a  /opt/intel/ipp/lib/intel64/libippvm.a)
endif()

if(FFTSCARF_FFT_FFTS)
    message("--   Enable FFTS (shipped with)")
    cat(${CMAKE_CURRENT_SOURCE_DIR}/include/fftscarf_ffts.h ${CMAKE_CURRENT_BINARY_DIR}/fftscarf.h)
    message("--       Configuring FFTS library...")
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_ffts.cpp)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ffts ${CMAKE_CURRENT_BINARY_DIR}/ffts_build)
#    add_custom_target(libffts
#                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ffts_build/
#                     COMMAND make)
#    add_dependencies(fftscarf libffts)
    if(WIN32)
        target_link_libraries(fftscarf ${CMAKE_CURRENT_BINARY_DIR}/ffts_build/Release/ffts_static.lib)
    else()
        target_link_libraries(fftscarf ${CMAKE_CURRENT_BINARY_DIR}/ffts_build/libffts.a)
    endif()
endif()

if(FFTSCARF_FFT_PFFFT)
    message("--   Enable Pretty Fast FFT (PFFFT) (built-in)")
    cat(${CMAKE_CURRENT_SOURCE_DIR}/include/fftscarf_pffft.h ${CMAKE_CURRENT_BINARY_DIR}/fftscarf.h)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_pffft.cpp)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/pffft/pffft.c)
endif()

if(FFTSCARF_FFT_FFTW3)
    message("--   Enable FFTW3 (external lib)")
    cat(${CMAKE_CURRENT_SOURCE_DIR}/include/fftscarf_fftw3.h ${CMAKE_CURRENT_BINARY_DIR}/fftscarf.h)
    include_directories(C:/projects/fftscarf/lib/libfft)
    if(FFTSCARF_PRECISION_SINGLE)
        if(WIN32)
            target_link_libraries(fftscarf C:/projects/fftscarf/lib/libfft/libfftw3f-3.lib)
        else()
            target_link_libraries(fftscarf fftw3f)
        endif()
    endif()
    if(FFTSCARF_PRECISION_DOUBLE)
        if(WIN32)
            target_link_libraries(fftscarf C:/projects/fftscarf/lib/libfft/libfftw3-3.lib)
        else()
            target_link_libraries(fftscarf fftw3)
        endif()
    endif()
    if(FFTSCARF_PRECISION_LONGDOUBLE)
        if(WIN32)
            target_link_libraries(fftscarf C:/projects/fftscarf/lib/libfft/libfftw3l-3.lib)
        else()
            target_link_libraries(fftscarf fftw3l)
        endif()
    endif()
endif()

if(FFTSCARF_FFT_OOURA)
    message("--   Enable Ooura (built-in)")
    cat(${CMAKE_CURRENT_SOURCE_DIR}/include/fftscarf_ooura.h ${CMAKE_CURRENT_BINARY_DIR}/fftscarf.h)
    if(FFTSCARF_PRECISION_DEFAULT EQUAL 32)
        add_definitions(-DOOFLOAT_SINGLE)
    endif()
    if(FFTSCARF_PRECISION_DEFAULT EQUAL 64)
        add_definitions(-DOOFLOAT_DOUBLE)
    endif()
    if(FFTSCARF_PRECISION_DEFAULT EQUAL 128)
        add_definitions(-DOOFLOAT_LONGDOUBLE)
    endif()
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_ooura.cpp)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ooura/fftsg.c)
#    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ooura/fft4g.cpp)
endif()

if(FFTSCARF_FFT_FFTREAL)
    message("--   Enable FFTReal (built-in)")
    cat(${CMAKE_CURRENT_SOURCE_DIR}/include/fftscarf_fftreal.h ${CMAKE_CURRENT_BINARY_DIR}/fftscarf.h)
endif()

# if(FFTSCARF_FFT_DJBFFT)
#     message("--   Enable DJBFFT (shipped with)")
# #    message("--    Configuring DJBFFT library")
#     target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_djbfft.cpp)
#     target_link_libraries(fftscarf ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/djbfft-0.76/djbfft.a)
# endif(FFTSCARF_FFT_DJBFFT)

if(FFTSCARF_FFT_DFT)
    message("--   Enable Standard DFT (built-in)")
    cat(${CMAKE_CURRENT_SOURCE_DIR}/include/fftscarf_dft.h ${CMAKE_CURRENT_BINARY_DIR}/fftscarf.h)
endif()


include_directories(${CMAKE_CURRENT_BINARY_DIR})
#include_directories(include)
include_directories(fftlibs)


# Testing ---------------------------------------------------------------------

option(FFTSCARF_BUILD_TEST "Build tests" OFF)
if(FFTSCARF_BUILD_TEST)
    message("\n-- FFTScarf configuring tests")
    message("-- Forcing FFTW3 for comparison purpose")
    add_subdirectory(test)
endif()

message("")


# Optimization -----------------------------------------------------------------

if(WIN32) # TODO Should be test on GCC because the followin might work on mingw
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Ox -D_USE_MATH_DEFINES /arch:SSE") # Needed by PFFFT
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox -D_USE_MATH_DEFINES /arch:SSE") # Needed by PFFFT
else() # TODO Should be test on GCC because the followin might work on mingw
    # Speed up compilation
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

    # -O2 -O3 -Ofast
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

# Install ---------------------------------------------------------------------
message("-- Installation path: " ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION include/
        FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include/
        FILES_MATCHING PATTERN "*.h")
install(TARGETS fftscarf DESTINATION lib)
