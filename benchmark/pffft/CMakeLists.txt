cmake_minimum_required (VERSION 2.8)

set(FFTNAME pffft)

set(EXEC benchmark_${FFTNAME})
project (${EXEC})
add_executable(${EXEC} benchmark.cpp)
target_link_libraries(${EXEC} ${Boost_LIBRARIES})
target_link_libraries(${EXEC} fftscarf)


add_custom_command(OUTPUT benchmark.log
                COMMAND ./${EXEC}
                DEPENDS ${EXEC}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(benchmark_run_${FFTNAME}
                    DEPENDS benchmark.log)

# Optimization -----------------------------------------------------------------

# Speed up compilation
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe")

# -O2 -O3 -Ofast
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast") # Almost no gain and not advised
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast") # Almost no gain and not advised
