cmake_minimum_required (VERSION 3.2)
project (fftscarf)

add_library(fftscarf src/fftscarf.cpp)

message("\n-- FFTScarf configuring lib")

option(FFTSCARF_PRECISION_SINGLE "Enable single precision" ON)
option(FFTSCARF_PRECISION_DOUBLE "Enable double precision" ON)
option(FFTSCARF_PRECISION_DEFAULTSINGLE "Set default precision to single" ON)

option(FFTSCARF_PLAN_PROTECTACCESS "Protect access of the plan resize and execution using a mutex (need boost library)." OFF)
if(FFTSCARF_PLAN_PROTECTACCESS)
    message("-- Protect plan access (resize and execution) using a mutex (using Boost library)")
    set(Boost_USE_STATIC_LIBS on)
    find_package(Boost COMPONENTS thread REQUIRED)
    add_definitions(-DBOOST_DATE_TIME_POSIX_TIME_STD_CONFIG)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries(fftscarf ${Boost_LIBRARIES})
else(FFTSCARF_PLAN_PROTECTACCESS)
    message("-- Plan access are not protected")
endif(FFTSCARF_PLAN_PROTECTACCESS)

if(FFTSCARF_PRECISION_SINGLE)
    message("-- Enable single precision (float 32b)")
else(FFTSCARF_PRECISION_SINGLE)
    set(FFTSCARF_PRECISION_DEFAULTSINGLE OFF CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_FFT_FFTS OFF CACHE BOOL "docstring" FORCE)
    set(FFTSCARF_FFT_PFFFT OFF CACHE BOOL "docstring" FORCE)
endif(FFTSCARF_PRECISION_SINGLE)
if(FFTSCARF_PRECISION_DOUBLE)
    message("-- Enable double precision (float 64b)")
else(FFTSCARF_PRECISION_DOUBLE)
    set(FFTSCARF_PRECISION_DEFAULTSINGLE ON CACHE BOOL "docstring" FORCE)
endif(FFTSCARF_PRECISION_DOUBLE)
if(FFTSCARF_PRECISION_DEFAULTSINGLE)
    message("-- Default is single precision")
else(FFTSCARF_PRECISION_DEFAULTSINGLE)
    message("-- Default is double precision")
endif(FFTSCARF_PRECISION_DEFAULTSINGLE)

message("-- FFT Implementation")
option(FFTSCARF_FFT_IPP "IPP" OFF)
option(FFTSCARF_FFT_FFTS "FFTS" OFF)
option(FFTSCARF_FFT_PFFFT "PFFFT" OFF)
option(FFTSCARF_FFT_FFTW3 "FFTW3" OFF)
option(FFTSCARF_FFT_OOURA "Ooura's FFT" ON)
option(FFTSCARF_FFT_FFTREAL "FFTReal" OFF)
#option(FFTSCARF_FFT_DJBFFT "DJBFFT" OFF)

if(FFTSCARF_FFT_IPP)
    message("--   Enable Intel Integrated Performance Primitives (Intel IPP) (external lib)")
    include_directories(/opt/intel/ipp/include)
    target_link_libraries(fftscarf /opt/intel/ipp/lib/intel64/libipps.a /opt/intel/ipp/lib/intel64/libippcore.a  /opt/intel/ipp/lib/intel64/libippvm.a)
endif(FFTSCARF_FFT_IPP)

if(FFTSCARF_FFT_FFTS)
    message("--   Enable FFTS (shipped with)")
    message("--       Configuring FFTS library...")
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_ffts.cpp)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ffts ${CMAKE_CURRENT_BINARY_DIR}/ffts_build)
    add_custom_target(libffts
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ffts_build/
                     COMMAND make)
    add_dependencies(fftscarf libffts)
    target_link_libraries(fftscarf ${CMAKE_CURRENT_BINARY_DIR}/ffts_build/libffts.a)
endif(FFTSCARF_FFT_FFTS)

if(FFTSCARF_FFT_PFFFT)
    message("--   Enable Pretty Fast FFT (PFFFT) (built-in)")
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_pffft.cpp)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/pffft/pffft.c)
endif(FFTSCARF_FFT_PFFFT)

if(FFTSCARF_FFT_FFTW3)
    message("--   Enable FFTW3 (external lib)")
    include_directories(C:/projects/fftscarf/lib/libfft)
    if(FFTSCARF_PRECISION_SINGLE)
        if(WIN32)
            target_link_libraries(fftscarf C:/projects/fftscarf/lib/libfft/libfftw3f-3.lib)
        else(WIN32)
            target_link_libraries(fftscarf fftw3f)
        endif(WIN32)
    endif(FFTSCARF_PRECISION_SINGLE)
    if(FFTSCARF_PRECISION_DOUBLE)
        if(WIN32)
            target_link_libraries(fftscarf C:/projects/fftscarf/lib/libfft/libfftw3-3.lib)
        else(WIN32)
            target_link_libraries(fftscarf fftw3)
        endif(WIN32)
    endif(FFTSCARF_PRECISION_DOUBLE)
endif(FFTSCARF_FFT_FFTW3)

if(FFTSCARF_FFT_OOURA)
    message("--   Enable Ooura (built-in)")
    if(FFTSCARF_PRECISION_DEFAULTSINGLE)
        add_definitions(-DOOFLOAT_FLOAT)
    endif(FFTSCARF_PRECISION_DEFAULTSINGLE)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_ooura.cpp)
    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ooura/fftsg.c)
#    target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/ooura/fft4g.cpp)
endif(FFTSCARF_FFT_OOURA)

if(FFTSCARF_FFT_FFTREAL)
    message("--   Enable FFTReal (built-in)")
endif(FFTSCARF_FFT_FFTREAL)

# if(FFTSCARF_FFT_DJBFFT)
#     message("--   Enable DJBFFT (shipped with)")
# #    message("--    Configuring DJBFFT library")
#     target_sources(fftscarf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/fftscarf_djbfft.cpp)
#     target_link_libraries(fftscarf ${CMAKE_CURRENT_SOURCE_DIR}/fftlibs/djbfft-0.76/djbfft.a)
# endif(FFTSCARF_FFT_DJBFFT)

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
                        COMMAND test_version)
    add_dependencies(test_run_version test_version)

    add_custom_target(test_run_fftlibs
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test
                     COMMAND test_fftlibs)
    add_dependencies(test_run_fftlibs test_fftlibs)

endif(FFTSCARF_BUILD_TEST)

message("")

# Optimization -----------------------------------------------------------------

if(WIN32) # TODO Should be test on GCC because the followin might work on mingw
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Ox -D_USE_MATH_DEFINES /arch:SSE") # Needed by PFFFT
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2")
else(WIN32) # TODO Should be test on GCC because the followin might work on mingw
    # Speed up compilation
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

    # -O2 -O3 -Ofast
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif(WIN32)

# Install ---------------------------------------------------------------------
message("-- Installation path: " ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ DESTINATION include/
        FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include/
        FILES_MATCHING PATTERN "*.h")
install(TARGETS fftscarf DESTINATION lib)
