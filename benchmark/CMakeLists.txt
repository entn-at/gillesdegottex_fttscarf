cmake_minimum_required (VERSION 3.2)
project (benchmark)

set(FFTSCARF_PRECISION_DEFAULT 32 CACHE STRING "docstring")
option(FFTSCARF_PRECISION_DEFAULT "Set default precision (32, 64, 128)")

option(FFTSCARF_FFT_IPP "" OFF)
option(FFTSCARF_FFT_FFTS "" ON)
option(FFTSCARF_FFT_PFFFT "" ON)
option(FFTSCARF_FFT_FFTW3 "" ON)
option(FFTSCARF_FFT_OOURA "" ON)
option(FFTSCARF_FFT_FFTREAL "" ON)
#option(FFTSCARF_FFT_DJBFFT "" OFF)
option(FFTSCARF_FFT_DFT "" OFF)

# Add Boost library
set(Boost_USE_STATIC_LIBS on)
find_package(Boost COMPONENTS chrono REQUIRED)
add_definitions(-DBOOST_DATE_TIME_POSIX_TIME_STD_CONFIG)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add FFTScarf library
set(FFTSCARF_BUILD_TEST OFF CACHE BOOL "docstring" FORCE)
set(FFTSCARF_PLAN_PROTECTACCESS OFF CACHE BOOL "docstring" FORCE)
if(FFTSCARF_PRECISION_DEFAULT EQUAL 32)
    message("\n-- Use single precision (float 32b)")
    set(FFTSCARF_PRECISION_SINGLE ON CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_PRECISION_DOUBLE OFF CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_PRECISION_LONGDOUBLE OFF CACHE BOOL "docstring" FORCE)
    add_definitions(-DOOFLOAT_SINGLE)
elseif(FFTSCARF_PRECISION_DEFAULT EQUAL 64)
    message("\n-- Use double precision (float 64b)")
    set(FFTSCARF_PRECISION_SINGLE OFF CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_PRECISION_DOUBLE ON CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_PRECISION_LONGDOUBLE OFF CACHE BOOL "docstring" FORCE)
    add_definitions(-DOOFLOAT_DOUBLE)
elseif(FFTSCARF_PRECISION_DEFAULT EQUAL 128)
    message("\n-- Use long double precision (float 128b)")
    set(FFTSCARF_PRECISION_SINGLE OFF CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_PRECISION_DOUBLE OFF CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_PRECISION_LONGDOUBLE ON CACHE BOOL "docstring" FORCE)
    add_definitions(-DOOFLOAT_LONGDOUBLE)
endif()
if(FFTSCARF_FFT_IPP)
    include_directories(/opt/intel/ipp/include)
endif()
include_directories(${CMAKE_SOURCE_DIR}/../fftlibs)
include_directories(${CMAKE_BINARY_DIR}/fftscarf)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_BINARY_DIR}/fftscarf)


# A make target for running the benchmarks
add_custom_target(benchmark_run DEPENDS benchmark_run_ooura)  # TODO rm benchmark_run_ooura

if(FFTSCARF_FFT_IPP)
    add_subdirectory(ipp)
    add_dependencies(benchmark_run benchmark_run_ipp)
endif(FFTSCARF_FFT_IPP)
if(FFTSCARF_FFT_FFTS)
    add_subdirectory(ffts)
    add_dependencies(benchmark_run benchmark_run_ffts)
endif(FFTSCARF_FFT_FFTS)
if(FFTSCARF_FFT_PFFFT)
    add_subdirectory(pffft)
    add_dependencies(benchmark_run benchmark_run_pffft)
endif(FFTSCARF_FFT_PFFFT)
if(FFTSCARF_FFT_FFTW3)
    add_subdirectory(fftw3)
    add_dependencies(benchmark_run benchmark_run_fftw3)
endif(FFTSCARF_FFT_FFTW3)
if(FFTSCARF_FFT_OOURA)
    add_subdirectory(ooura)
    add_dependencies(benchmark_run benchmark_run_ooura)
endif(FFTSCARF_FFT_OOURA)
if(FFTSCARF_FFT_FFTREAL)
    add_subdirectory(fftreal)
    add_dependencies(benchmark_run benchmark_run_fftreal)
endif(FFTSCARF_FFT_FFTREAL)
if(FFTSCARF_FFT_DFT)
    add_subdirectory(dft)
    add_dependencies(benchmark_run benchmark_run_dft)
endif(FFTSCARF_FFT_DFT)

# A make target for ploting the results of the benchmark
add_custom_target(plot
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
                    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/plot.py)
