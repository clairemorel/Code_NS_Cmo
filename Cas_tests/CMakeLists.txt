#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

set(EXECUTABLE_OUTPUT_PATH ../)

include_directories(../include)

link_directories(/lib/${CMAKE_BUILD_TYPE})

file(GLOB CODE_src1 main_NS.cpp)
file(GLOB CODE_src2 main_chaleur_stat.cpp)
file(GLOB CODE_src3 main_chaleur_cart.cpp)

add_executable(exe_NS ${CODE_src1})
add_executable(exe_chaleur_stat ${CODE_src2})
add_executable(exe_chaleur_cart ${CODE_src3})

target_link_libraries(exe_NS NavierStokes)
target_link_libraries(exe_chaleur_stat NavierStokes)
target_link_libraries(exe_chaleur_cart NavierStokes)

target_link_libraries(exe_NS ${BITPIT_LIBRARIES})
target_link_libraries(exe_chaleur_stat ${BITPIT_LIBRARIES})
target_link_libraries(exe_chaleur_cart ${BITPIT_LIBRARIES})

target_link_libraries(exe_NS ${PETSC_LIBRARIES})
target_link_libraries(exe_chaleur_stat ${PETSC_LIBRARIES})
target_link_libraries(exe_chaleur_cart ${PETSC_LIBRARIES})