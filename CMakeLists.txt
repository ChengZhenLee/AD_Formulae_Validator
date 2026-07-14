cmake_minimum_required(VERSION 3.14)
project(ADValidator)

# Use C++ 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the executable with all source files
option(BUILD_TESTS "Build unit tests" ON)

# GoogleTest
if(BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
    )
    FetchContent_MakeAvailable(googletest)
endif()

# Build a generator library from the generator folder.
file(GLOB_RECURSE GENERATOR_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/generator/*.cpp"
)
# Exclude adDrivers.cpp/adHelper.cpp
list(FILTER GENERATOR_SRCS EXCLUDE REGEX "/(adDrivers|adHelper)\\.cpp$")
add_library(generator_lib STATIC ${GENERATOR_SRCS})

# Copy the generator and include folders into the build directory.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/generator DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(generator_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/generator
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ad
)

# Main executable: main.cpp is just argument parsing + stage orchestration;
# the pipeline stages and all terminal-output formatting live in cli/.
add_executable(ADValidator
    main.cpp
    cli/console.cpp
    cli/pipeline.cpp
)
target_link_libraries(ADValidator PRIVATE generator_lib)

target_include_directories(ADValidator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Assemble a self-contained redistributable folder: the executable plus the
# resource files it needs.
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

# Unit tests
if(BUILD_TESTS)
    enable_testing()

    file(GLOB_RECURSE TEST_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
    add_executable(generator_tests ${TEST_SRCS})
    target_link_libraries(generator_tests PRIVATE generator_lib GTest::gtest_main)
    target_include_directories(generator_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    include(GoogleTest)
    gtest_discover_tests(generator_tests)
endif()