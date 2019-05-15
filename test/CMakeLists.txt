macro(buildTarget TARGET SRC ADD_TO_TESTSUITE)
    message(STATUS "Test lib: " ${TEST_BUILD_LIBRARY})
    add_executable(${TARGET} ${SRC})
    target_link_libraries(${TARGET} ${TEST_BUILD_LIBRARY} ${GTEST_LIBRARIES})
    if (${ADD_TO_TESTSUITE})
        add_test(${TARGET} ${TARGET})
    endif()
endmacro(buildTarget)

buildTarget(testGraph TestGraph.cpp ON)
buildTarget(testDynamicBitset TestDynamicBitset.cpp ON)

