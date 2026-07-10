cmake_minimum_required(VERSION 3.10)
project(ADValidator)

# Use C++ 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the executable with all source files
option(BUILD_TESTS "Build unit tests" ON)

# Build a generator library from the generator folder. Recursive because
# generator/ groups its source into subfolders by role (generator/codegen,
# generator/symbolic) - see the "Project layout" section of the README.
file(GLOB_RECURSE GENERATOR_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/generator/*.cpp"
)
# Exclude adDrivers.cpp/adHelper.cpp: not source, but ADValidator's own
# runtime-generated output (each has its own main()), left sitting flat in
# generator/ from a previous run. Without this, re-running cmake after
# running ADValidator once would sweep them into this library too.
list(FILTER GENERATOR_SRCS EXCLUDE REGEX "/(adDrivers|adHelper)\\.cpp$")
add_library(generator_lib STATIC ${GENERATOR_SRCS})

# Copy the generator and include folders into the build directory. At
# runtime ADValidator resolves both as siblings of its own executable (see
# getExecutableDir() in main.cpp), so this mirrors the layout it expects.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/generator DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(generator_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/generator
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ad
)

# Main executable: main.cpp is just argument parsing + stage orchestration;
# the pipeline stages and all terminal-output formatting live in cli/ (see
# cli/pipeline.* and cli/console.*).
add_executable(ADValidator
    main.cpp
    cli/console.cpp
    cli/pipeline.cpp
)
target_link_libraries(ADValidator PRIVATE generator_lib)

# So `#include "generator/..."` and `#include "cli/..."` resolve the same
# way from every translation unit above, regardless of which subdirectory
# it's compiled from.
target_include_directories(ADValidator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Assemble a self-contained redistributable folder: the executable plus the
# resource files it needs (and its full generator/ source, for a peer who
# wants to read it, not just run it). Hand peers this dist/ folder instead
# of the whole build tree.
#
# generator/ is mirrored wholesale with copy_directory (like include/
# already is) rather than an explicit file list, so this never needs
# hand-updating again as source moves between generator/'s subfolders. The
# explicit removals right after strip out ADValidator's own runtime output
# (adDrivers.*/adHelper.*/*.bin/logs/etc - see .gitignore for the full set)
# in case a previous run left any sitting in the source tree being copied
# from; -f suppresses the "file doesn't exist" error since normally none of
# them do.
set(DIST_DIR ${CMAKE_CURRENT_BINARY_DIR}/dist)
add_custom_command(TARGET ADValidator POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DIST_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:ADValidator> ${DIST_DIR}/
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/generator ${DIST_DIR}/generator
    COMMAND ${CMAKE_COMMAND} -E remove -f
        ${DIST_DIR}/generator/adDrivers.h
        ${DIST_DIR}/generator/adDrivers.cpp
        ${DIST_DIR}/generator/adDrivers${CMAKE_EXECUTABLE_SUFFIX}
        ${DIST_DIR}/generator/adDrivers_build.log
        ${DIST_DIR}/generator/adHelper.h
        ${DIST_DIR}/generator/adHelper.cpp
        ${DIST_DIR}/generator/adHelper${CMAKE_EXECUTABLE_SUFFIX}
        ${DIST_DIR}/generator/adHelper_build.log
        ${DIST_DIR}/generator/parameters.bin
        ${DIST_DIR}/generator/derivatives_seeds.bin
        ${DIST_DIR}/generator/derivatives.bin
        ${DIST_DIR}/generator/results.bin
        ${DIST_DIR}/generator/equations.txt
        ${DIST_DIR}/generator/validation_results.txt
        ${DIST_DIR}/generator/.pipeline_cache
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/include/ad ${DIST_DIR}/include/ad
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/include/nlohmann ${DIST_DIR}/include/nlohmann
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/include/Eigen ${DIST_DIR}/include/Eigen
    COMMENT "Assembling redistributable package in ${DIST_DIR}"
)