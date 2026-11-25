cmake_minimum_required(VERSION 3.16)
project(computational-geometry-cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(geometry
  src/closest_pair.cpp
  src/neighbors.cpp
  src/convex_hull.cpp
  src/delaunay.cpp
  src/shoelace.cpp
  src/segments.cpp
  src/bounding_box.cpp
  src/laplacian.cpp
  src/svg.cpp
)
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(geometry PRIVATE -Wall -Wextra)

add_executable(geometry-cli src/cli.cpp)
target_link_libraries(geometry-cli PRIVATE geometry)
set_target_properties(geometry-cli PROPERTIES OUTPUT_NAME geometry)
target_compile_options(geometry-cli PRIVATE -Wall -Wextra)

add_executable(geometry_tests tests/test_geometry.cpp)
target_link_libraries(geometry_tests PRIVATE geometry)
target_compile_options(geometry_tests PRIVATE -Wall -Wextra)

enable_testing()
add_test(NAME geometry COMMAND geometry_tests)
