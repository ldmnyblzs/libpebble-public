cmake_minimum_required(VERSION 3.25)
project(libpebble)

include(CTest)
set(CMAKE_CXX_STANDARD 20)

include(${CMAKE_BINARY_DIR}/conan_paths.cmake OPTIONAL)
find_package(Boost 1.81.0 REQUIRED)
find_package(CGAL 5.5.1 REQUIRED)
find_package(bliss 0.77 REQUIRED)
find_package(Eigen3 3.3.9 REQUIRED)

add_library(pebble INTERFACE)
target_compile_definitions(pebble INTERFACE BOOST_ALLOW_DEPRECATED_HEADERS=1)
target_include_directories(pebble INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(pebble SYSTEM INTERFACE ${Boost_INCLUDE_DIRS} ${bliss_INCLUDE_DIRS})
target_link_libraries(pebble INTERFACE CGAL::CGAL ${bliss_LIBRARIES} Eigen3::Eigen)

add_subdirectory(example)
add_subdirectory(test)

find_package(Doxygen REQUIRED)
set(DOXYGEN_EXCLUDE_SYMBOLS "pebble::internal")
set(DOXYGEN_EXCLUDE_PATTERNS "utility.hpp")
doxygen_add_docs(docs include/pebble example ALL)
