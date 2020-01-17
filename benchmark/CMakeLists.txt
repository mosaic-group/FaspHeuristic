project(FaspTightCutBenchmark VERSION 1.0.0 LANGUAGES C CXX)
message(STATUS "FASP: Building benchmarks")

find_package(HDF5 REQUIRED)

add_executable(benchmark main.cpp tools/easylogging++.cc)

target_link_libraries(benchmark PRIVATE ${HDF5_LIBRARIES})
target_compile_features(benchmark PRIVATE cxx_std_17)
if(CMAKE_COMPILER_IS_GNUCC)
    target_compile_options(benchmark PUBLIC $<$<CONFIG:Debug>:-O0 -g>)
    target_compile_options(benchmark PUBLIC $<$<CONFIG:Release>:-O4 -pthread -DNDEBUG>)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(benchmark PUBLIC $<$<CONFIG:Debug>:-O0 -g>)
    target_compile_options(benchmark PUBLIC $<$<CONFIG:Release>:-Ofast -DNDEBUG>)
endif()
