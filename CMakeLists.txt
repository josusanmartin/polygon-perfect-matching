cmake_minimum_required(VERSION 3.16)
project(mwpm CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

include(FetchContent)
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(json)

add_library(mwpm_lib
    src/generator.cpp
    src/geometry.cpp
    src/io.cpp
    src/algorithms/brute_force.cpp
    src/algorithms/dynamic_programming.cpp
    src/algorithms/marcotte_suri.cpp
)

target_include_directories(mwpm_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(mwpm_lib PUBLIC nlohmann_json::nlohmann_json)

add_executable(mwpm src/main.cpp)
target_link_libraries(mwpm PRIVATE mwpm_lib)
