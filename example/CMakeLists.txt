add_executable(pebble-cli pebble-cli.cpp)
target_link_libraries(pebble-cli PRIVATE pebble)

add_executable(graph2tikz graph2tikz.cpp)
target_link_libraries(graph2tikz PRIVATE pebble)
