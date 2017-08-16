project(boost_graph_tests)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "-Wall -g -std=c++14")

find_package(Boost 1.64.0 REQUIRED)
find_package(OpenCV REQUIRED)

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(boost_graph_tests ${OpenCV_LIBS})
