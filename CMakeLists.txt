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

# Copy the generator folder into the build
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/generator DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(generator_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/generator
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ad
)

# Main executable (minimal runner for generator-focused builds)
add_executable(ADValidator main.cpp)
target_link_libraries(ADValidator PRIVATE generator_lib)