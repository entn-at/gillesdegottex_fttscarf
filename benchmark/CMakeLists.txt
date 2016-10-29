cmake_minimum_required (VERSION 2.8)
project (eval_fft)

option(FFTSCARF_PRECISION_FLOAT32 "" ON)
option(FFTSCARF_BUILD_TEST "" OFF)
option(FFTSCARF_FFT_OOURA "" ON)
option(FFTSCARF_FFT_FFTREAL "" ON)
option(FFTSCARF_FFT_PFFFT "" OFF)
option(FFTSCARF_FFT_DJBFFT "" OFF)
option(FFTSCARF_FFT_FFTS "" OFF)
option(FFTSCARF_FFT_FFTW3 "" OFF)
option(FFTSCARF_FFT_IPP "" OFF)

# Add Boost library
find_package(Boost COMPONENTS system filesystem chrono REQUIRED)
add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_DATE_TIME_POSIX_TIME_STD_CONFIG)
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add Voimoo library
if(FFTSCARF_PRECISION_FLOAT32)
    add_definitions(-DOOFLOAT_FLOAT) # TODO Annoying to have to redefine this here
endif(FFTSCARF_PRECISION_FLOAT32)
include_directories(${CMAKE_SOURCE_DIR}/../include/)
include_directories(${CMAKE_BINARY_DIR}/fftscarf/include)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_BINARY_DIR}/fftscarf)

if(FFTSCARF_FFT_OOURA)
    add_subdirectory(ooura)
endif(FFTSCARF_FFT_OOURA)
if(FFTSCARF_FFT_FFTREAL)
    add_subdirectory(fftreal)
endif(FFTSCARF_FFT_FFTREAL)
if(FFTSCARF_FFT_PFFFT)
    add_subdirectory(pffft)
endif(FFTSCARF_FFT_PFFFT)
if(FFTSCARF_FFT_DJBFFT)
    add_subdirectory(djbfft)
endif(FFTSCARF_FFT_DJBFFT)
if(FFTSCARF_FFT_FFTS)
    add_subdirectory(ffts)
endif(FFTSCARF_FFT_FFTS)
if(FFTSCARF_FFT_FFTW3)
    add_subdirectory(fftw3)
endif(FFTSCARF_FFT_FFTW3)
if(FFTSCARF_FFT_IPP)
    add_subdirectory(ipp)
    include_directories(/opt/intel/ipp/include)
endif(FFTSCARF_FFT_IPP)

# A make target for running the benchmarks
add_custom_target(benchmark_run
                    DEPENDS benchmark_run_ooura benchmark_run_fftreal) # benchmark_run_fftw3 benchmark_run_ffts benchmark_run_djbfft benchmark_run_pffft benchmark_run_ipp

# A make target for ploting the results of the benchmark
add_custom_target(plot
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
                    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/plot.py)
