###############################################################################
# FASP heuristic
###############################################################################
cmake_minimum_required(VERSION 3.10)
project(FaspTightCut VERSION 1.0.0 LANGUAGES C CXX)

#project options
option(FASP_HEURISTIC_BUILD_TESTS "Build tests?" ON)
option(FASP_HEURISTIC_BUILD_BENCHMARKS "Build benchmarks?" OFF)
option(FASP_HEURISTIC_INSTALL "Install FaspHeuristic?" OFF)

message(STATUS "FASP build type: " ${CMAKE_BUILD_TYPE})
###############################################################################
# Create FASP TightCut library
###############################################################################
set(LIB_NAMESPACE FaspTightCut)
set(LIB_NAME FaspTightCut)

add_library(${LIB_NAME} INTERFACE)
target_include_directories(${LIB_NAME} INTERFACE
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

###############################################################################
# Install
###############################################################################
message(STATUS "FASP: Install library in [${CMAKE_INSTALL_PREFIX}]")

include(GNUInstallDirs)

# copy includes
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# generate target target file
set(TARGET_FILE_NAME ${LIB_NAME}Targets)
install(TARGETS ${LIB_NAME} EXPORT ${TARGET_FILE_NAME} INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT ${LIB_NAME}Targets  NAMESPACE ${LIB_NAMESPACE}:: DESTINATION lib/cmake/${LIB_NAME})

# generate version config version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(${LIB_NAME}ConfigVersion.cmake
        VERSION ${CMAKE_PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
        )

# generate config file
configure_file(cmake/${LIB_NAME}Config.cmake.in ${LIB_NAME}Config.cmake @ONLY)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake
        DESTINATION
        lib/cmake/${LIB_NAME})

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(${LIB_NAMESPACE}::${LIB_NAME} ALIAS ${LIB_NAME})


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
# Lib needed for tests and benchmarks
###############################################################################
add_library(testObjLib OBJECT src/tools/easylogging++.cc)
set(FASP_BENCHMARK_LIB benchmarkLib)
add_library(${FASP_BENCHMARK_LIB} STATIC  $<TARGET_OBJECTS:testObjLib> )
target_link_libraries(${FASP_BENCHMARK_LIB} ${LIB_NAME})

###############################################################################
# Benchmarks
###############################################################################
if (FASP_HEURISTIC_BUILD_BENCHMARKS)
    message(STATUS "FASP: Building benchmarks")
    add_executable(benchmark src/main.cpp)
    target_link_libraries(benchmark PRIVATE ${FASP_BENCHMARK_LIB} ${HDF5_LIBRARIES})
endif()


###############################################################################
# Tests
###############################################################################
if (FASP_HEURISTIC_BUILD_TESTS)
    message(STATUS "FASP: Building tests")
    set(TEST_BUILD_LIBRARY ${FASP_BENCHMARK_LIB})
    enable_testing()
    add_subdirectory(test)
endif()


#message(STATUS "---------------------- VARS BEG ---------------------")
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#message(STATUS "---------------------- VARS END ---------------------")
