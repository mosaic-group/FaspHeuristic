cmake_minimum_required(VERSION 3.10)
project(FaspHeuristic)

set(CMAKE_CXX_STANDARD 17)

# Find needed package
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

# Should link time optimization be used? (some compilers have problems with it)
set(USE_LTO OFF)

# Setup all needed compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -pg -std=c++17")
if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS_RELEASE "-O4 -DNDEBUG -pthread -g")
    if(${USE_LTO})
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -pthread")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
    if(${USE_LTO})
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
    endif()

    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
endif()

# Define sources
include_directories(src)
include_directories(src/tools)
set(SOURCE_FILES
        src/tools/easylogging++.cc
)

# And build output
add_library(testObjLib OBJECT ${SOURCE_FILES})
set(STATIC_TARGET_NAME staticLib)

add_library(${STATIC_TARGET_NAME} STATIC $<TARGET_OBJECTS:testObjLib>)
set(TEST_BUILD_LIBRARY ${STATIC_TARGET_NAME})

add_executable(benchmark src/main.cpp)
target_link_libraries(benchmark PRIVATE ${TEST_BUILD_LIBRARY} ${HDF5_LIBRARIES})

add_executable(playground src/playground.cpp)
target_link_libraries(playground PRIVATE ${TEST_BUILD_LIBRARY} ${HDF5_LIBRARIES})


# Build tests
message(STATUS "Building tests")
set(BUILD_GMOCK ON CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
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
