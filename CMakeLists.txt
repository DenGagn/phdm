# --------------------------------------------------------------------------- #
# Author:       Denis Gagnon                  <denisg6@hotmail.com>           #
# Date created: 2018-05-01                                                    #
# Description:  CMake compilation instructions for phdm                       #
# ----------------------------------------------------------------------------#

project(phdm)

# -- CMake version and installation directory.
# CMake version
cmake_minimum_required(VERSION 3.1)

# C++14 Standard required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX /usr)
endif()
LIST (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})

# -- Required dependency: Boost.
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
#SET (LIBS ${LIBS} ${Boost_LIBRARIES} )

# -- Required dependency: Armadillo
find_package(armadillo REQUIRED)
include_directories(${armadillo_INCLUDE_DIRS})
set(LIBS ${LIBS} ${armadillo_LIBRARIES})

# -- Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# -- Compiler config for RELEASE build.
# Configuration for the GCC compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3")
  set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}    -O3")

elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3")
  set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}    -O3")
endif()

# -- Compiler config for DEBUG build.
# Configuration for the GCC compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -pg -g -Wall")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 -pg -g -Wall")

# Configuration for the Intel compiler.
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -debug all")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 -g -debug all")
endif()

# ----------------------------------------------------------------- #
# --                  Compilation Instructions                   -- #
# ----------------------------------------------------------------- #
# -- Included files
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)

# -- Install mellotron headers
install (DIRECTORY include/          DESTINATION include)

# -- Output binaries in directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/simulations)

# -- Rules for every target

add_executable("GrapheneSine.o" "simulations/GrapheneSine.cpp")
target_link_libraries("GrapheneSine.o" ${LIBS})
install(TARGETS "GrapheneSine.o" RUNTIME DESTINATION bin)

add_executable("GrapheneExp.o" "simulations/GrapheneExp.cpp")
target_link_libraries("GrapheneExp.o" ${LIBS})
install(TARGETS "GrapheneExp.o" RUNTIME DESTINATION bin)

