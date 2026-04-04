cmake_minimum_required(VERSION 3.10)
project(ADValidator)

# Use C++ 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the executable with all source files
add_executable(ADValidator
    main.cpp
    lighthouse/1/fd.cpp
    lighthouse/1/va.cpp
    lighthouse/1/vt.cpp
    lighthouse/2/vasa.cpp
    lighthouse/2/vast.cpp
    lighthouse/2/vava.cpp
    lighthouse/2/vavt.cpp
    lighthouse/2/vtva.cpp
    lighthouse/2/vtvt.cpp
    lighthouse/3/vasasa.cpp
    lighthouse/3/vasast.cpp
    lighthouse/3/vastsa.cpp
    lighthouse/3/vastst.cpp
    lighthouse/3/vtvasa.cpp
    lighthouse/3/vtvast.cpp
    lighthouse/3/vtvtva.cpp
    lighthouse/3/vtvtvt.cpp
    lighthouse/4/vtvtvtvt.cpp
)

# Include the headers
target_include_directories(ADValidator PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ad
    ${CMAKE_CURRENT_SOURCE_DIR}/lighthouse
)