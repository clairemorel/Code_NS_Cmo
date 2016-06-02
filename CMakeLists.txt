#------------------------------------------------------------------------------------#
# CMake settings
#------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6)

project("Code_parallele")

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#

set(BITPIT_DIR "/home/cmorel/bitpit-private/lib/cmake/bitpit-1.2" CACHE PATH "The directory containig the bitpit library")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

set(PETSC_DIR $ENV{PETSC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})
set(CMAKE_CXX_FLAGS "-g -O0")

#------------------------------------------------------------------------------------#
# External dependencies
#------------------------------------------------------------------------------------#

find_package(MPI)

find_package(BITPIT REQUIRED)

include(FindPETSc)

# If BITPIT was not found, explain to the user how to specify its location.
if (NOT BITPIT_FOUND)
    set(BITPIT_DIR_MESSAGE "BITPIT not found. Set the BITPIT_DIR cmake cache entry to the directory containing BITPITConfig.cmake")

    if (BITPIT_FIND_REQUIRED)
        message(FATAL_ERROR ${BITPIT_DIR_MESSAGE})
    elseif (NOT BITPIT_FIND_QUIETLY)
        message(STATUS ${BITPIT_DIR_MESSAGE})
    endif ()
endif ()

add_definitions(-std=c++0x)
add_definitions(${PETSC_DEFINITIONS})

include_directories(include)
include_directories(${PETSC_INCLUDES})

include(${BITPIT_USE_FILE})

add_subdirectory(Cas_tests)

file (GLOB_RECURSE source_files source/* include/*)

add_library(NavierStokes SHARED ${source_files})

# message("Bitpit_include_dir = ${BITPIT_INCLUDE_DIRS}")
# message("Bitpit_library_dir = ${BITPIT_LIBRARY_DIRS}")
# message("Bitpit_definitions = ${BITPIT_DEFINITIONS}")
# message("BITPIT_libraries=${BITPIT_LIBRARIES}")
# message("BITPIT_use= ${BITPIT_USE_FILE}")

# message("Petsc_dir = ${PETSC_DIR}")
# message("trouvé petsc ? = ${PETSC_FOUND}")
# message("Petsc_include = ${PETSC_INCLUDES} ")
# message("Petsc_definitions = ${PETSC_DEFINITIONS} ")
