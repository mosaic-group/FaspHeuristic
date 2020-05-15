project(FaspTightCutTools VERSION 1.0.0 LANGUAGES C CXX)

message(STATUS "FASP: Building tools")

# =============================================================================
# Build tools
# =============================================================================

macro(buildTarget TARGET SRC)
    add_executable(${TARGET} ${SRC})
    target_link_libraries(${TARGET} ${LIB_NAMESPACE}::${LIB_NAME})
    if (FASP_HEURISTIC_INSTALL)
        install(TARGETS ${TARGET} RUNTIME DESTINATION bin)
    endif()
endmacro(buildTarget)

buildTarget(tightCut tightCut.cpp)
buildTarget(tightCutWeighted tightCutWeighted.cpp)
