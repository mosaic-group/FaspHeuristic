project(FaspTightCutTest VERSION 1.0.0 LANGUAGES C CXX)

message(STATUS "FASP: Building tests")

# =============================================================================
# Download and build googletest
# =============================================================================

set(GOOGLETEST_DIRECTORY ${PROJECT_BINARY_DIR}/googletest)

configure_file(cmake/GoogleTestCMakeLists.txt.in ${GOOGLETEST_DIRECTORY}/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . RESULT_VARIABLE result WORKING_DIRECTORY ${GOOGLETEST_DIRECTORY})
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build . RESULT_VARIABLE result WORKING_DIRECTORY ${GOOGLETEST_DIRECTORY})
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Add googletest directly to build. This defines the gtest and gtest_main targets.
add_subdirectory(${GOOGLETEST_DIRECTORY}/src ${GOOGLETEST_DIRECTORY}/build EXCLUDE_FROM_ALL)


# =============================================================================
# Build tests
# =============================================================================

macro(buildTarget TARGET SRC ADD_TO_TESTSUITE)
    add_executable(${TARGET} ${SRC})
    target_link_libraries(${TARGET}  gtest_main gmock_main ${LIB_NAMESPACE}::${LIB_NAME})
    if (${ADD_TO_TESTSUITE})
        add_test(${TARGET} ${TARGET})
    endif()
endmacro(buildTarget)

buildTarget(testGraph TestGraph.cpp ON)
buildTarget(testDynamicBitset TestDynamicBitset.cpp ON)
buildTarget(testGraphIO TestGraphIO.cpp ON)
buildTarget(testStack TestStack.cpp ON)
buildTarget(testGraphFasp TestGraphFasp.cpp ON)
