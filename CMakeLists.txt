cmake_minimum_required(VERSION 3.10)
project(ADValidator)

# Use C++ 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the executable with all source files
option(BUILD_TESTS "Build unit tests" ON)

# Build a generator library from the generator folder
file(GLOB GENERATOR_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/generator/*.cpp"
)
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

# Main executable (minimal runner for generator-focused builds)
add_executable(ADValidator main.cpp)
target_link_libraries(ADValidator PRIVATE generator_lib)

# Assemble a minimal, self-contained redistributable folder: the executable
# plus only the resource files it actually needs at runtime (headers used by
# the code it generates and compiles, not the generator's own build
# artifacts/binaries/logs). Hand peers this dist/ folder instead of the
# whole build tree.
set(DIST_DIR ${CMAKE_CURRENT_BINARY_DIR}/dist)
add_custom_command(TARGET ADValidator POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DIST_DIR}/generator
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:ADValidator> ${DIST_DIR}/
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/generator/structures.h
        ${CMAKE_CURRENT_SOURCE_DIR}/generator/utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/generator/utils.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/generator/generator.h
        ${CMAKE_CURRENT_SOURCE_DIR}/generator/user_function.h
        ${CMAKE_CURRENT_SOURCE_DIR}/generator/readWrite.h
        ${CMAKE_CURRENT_SOURCE_DIR}/generator/configManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/generator/configs.txt
        ${DIST_DIR}/generator/
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/include/ad ${DIST_DIR}/include/ad
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/include/nlohmann ${DIST_DIR}/include/nlohmann
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/include/Eigen ${DIST_DIR}/include/Eigen
    COMMENT "Assembling redistributable package in ${DIST_DIR}"
)