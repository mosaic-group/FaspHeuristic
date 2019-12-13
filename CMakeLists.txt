###############################################################################
# FASP heuristic
###############################################################################
cmake_minimum_required(VERSION 3.10)
project(FaspHeuristic)

#project options
option(FASP_HEURISTIC_BUILD_TESTS "Build tests?" ON)
option(FASP_HEURISTIC_BUILD_BENCHMARKS "Build benchmarks?" ON)
option(FASP_HEURISTIC_INSTALL "Install FaspHeuristic?" OFF)

###############################################################################
# Find all required libraries
###############################################################################
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

###############################################################################
# Configure compiler options
###############################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS_RELEASE "-O4 -DNDEBUG -pthread -g")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -pthread")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
endif()

# Include paths
include_directories(src)
include_directories(src/tools)


###############################################################################
# Benchmarks
###############################################################################
add_library(testObjLib OBJECT src/tools/easylogging++.cc)
set(FASP_BENCHMARK_LIB benchmarkLib)
add_library(${FASP_BENCHMARK_LIB} STATIC $<TARGET_OBJECTS:testObjLib>)
add_executable(benchmark src/main.cpp)
target_link_libraries(benchmark PRIVATE ${FASP_BENCHMARK_LIB} ${HDF5_LIBRARIES})


###############################################################################
# Tests
###############################################################################
message(STATUS "Building tests")
set(BUILD_GMOCK ON CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(TEST_BUILD_LIBRARY ${FASP_BENCHMARK_LIB})
add_subdirectory("libs/googletest")
include_directories(libs/googletest/googlemock/include)
set(GTEST_LIBRARIES gtest gmock_main)
enable_testing()
add_subdirectory(test)


#message(STATUS "---------------------- VARS BEG ---------------------")
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#message(STATUS "---------------------- VARS END ---------------------")
