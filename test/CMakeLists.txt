find_package(GTest 1.12.1 REQUIRED)

add_executable(face_center_test face_center.cpp)
target_link_libraries(face_center_test PUBLIC GTest::gtest GTest::gtest_main PRIVATE pebble)

add_test(NAME face_center_test COMMAND face_center_test)
