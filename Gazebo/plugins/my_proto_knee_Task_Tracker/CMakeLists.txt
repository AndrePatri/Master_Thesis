cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(proto_controller)

find_package(gazebo REQUIRED)

include_directories( ${GAZEBO_INCLUDE_DIRS})
include_directories(blasfeo/include)
include_directories(hpipm/include)
add_subdirectory(blasfeo)
add_subdirectory(hpipm)
file(GLOB srcs "hpipm/include/*.h" "blasfeo/include/*.h")

link_directories(${GAZEBO_LIBRARY_DIRS})
	
	     
list(APPEND CMAKE_CXX_FLAGS "${GAZEBO_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GAZEBO_CXX_FLAGS}")

add_library(proto_controller SHARED proto_controller.cc)


target_link_libraries(proto_controller ${GAZEBO_LIBRARIES} blasfeo hpipm)

