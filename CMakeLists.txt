###############################################################################
# FASP heuristic
###############################################################################
cmake_minimum_required(VERSION 3.10)
project(FaspTightCut VERSION 1.0.0 LANGUAGES C CXX)

#project options
option(FASP_HEURISTIC_INSTALL "Install FaspHeuristic?" ON)
option(FASP_HEURISTIC_BUILD_TOOLS "Build tools?" ON)
option(FASP_HEURISTIC_BUILD_TESTS "Build tests?" OFF)
option(FASP_HEURISTIC_BUILD_BENCHMARKS "Build benchmarks?" OFF)


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
target_compile_options(${LIB_NAME} INTERFACE -Werror -Wall -pedantic)
target_compile_features(${LIB_NAME} INTERFACE cxx_std_17)

#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(${LIB_NAMESPACE}::${LIB_NAME} ALIAS ${LIB_NAME})

###############################################################################
# Install FASP TightCut
###############################################################################
if (FASP_HEURISTIC_INSTALL)
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
endif()

###############################################################################
# Tools
###############################################################################
if (FASP_HEURISTIC_BUILD_TOOLS)
    add_subdirectory(tools)
endif()

###############################################################################
# Benchmarks
###############################################################################
if (FASP_HEURISTIC_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()

###############################################################################
# Tests
###############################################################################
if (FASP_HEURISTIC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
