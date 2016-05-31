#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

set(EXECUTABLE_OUTPUT_PATH ../)

include_directories(../include)

link_directories(/lib/${CMAKE_BUILD_TYPE})

file(GLOB CODE_src main_NS.cpp)

add_executable(exe ${CODE_src})

target_link_libraries(exe NavierStokes)

target_link_libraries(exe ${BITPIT_LIBRARIES})

target_link_libraries(exe ${PETSC_LIBRARIES})