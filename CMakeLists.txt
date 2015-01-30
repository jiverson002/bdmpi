cmake_minimum_required(VERSION 2.8)
project(BDMP)

set(GKLIB_PATH "GKlib" CACHE PATH "path to GKlib")
set(SHARED FALSE CACHE BOOL "build a shared library")

set(HINSTALL_PATH "bdmpi" CACHE PATH "header installation subdirectory")
set(LINSTALL_PATH "bdmpi" CACHE PATH "library installation subdirectory")

if(MSVC)
  set(BDMP_INSTALL FALSE)
else()
  set(BDMP_INSTALL TRUE)
endif()

# Configure library type.
if(SHARED)
  set(BDMP_LIBRARY_TYPE SHARED)
else()
  set(BDMP_LIBRARY_TYPE STATIC)
endif(SHARED)

include(${GKLIB_PATH}/GKlibSystem.cmake)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -p")

# Collect the source files in the common directory
set(COMMON_PATH "common" CACHE PATH "path to common")
file(GLOB common_sources ${COMMON_PATH}/*.c)

# Add include directories.
include_directories(${GKLIB_PATH})
include_directories(include)

# Recursively look for CMakeLists.txt in subdirs.
add_subdirectory("GKlib")
add_subdirectory("include")
add_subdirectory("libbdmp")
add_subdirectory("bdmprun")
add_subdirectory("test")
add_subdirectory("scripts")

add_subdirectory("apps/bdmp")           # DISTRM  this is for my grep -v script
add_subdirectory("apps/bdmp-oc")        # DISTRM
add_subdirectory("apps/oc")             # DISTRM
add_subdirectory("apps/mpi")            # DISTRM
add_subdirectory("apps/misc")           # DISTRM
