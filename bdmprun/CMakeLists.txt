# These programs use internal bdmp data structures.
#include_directories(../libbdmp)

# Find sources.
file(GLOB bdmp_sources *.c)

# Build program.
add_executable(bdmprun ${common_sources} ${bdmp_sources})
#add_dependencies(bdmprun ${GKlib_includes})
foreach(prog bdmprun)
  target_link_libraries(${prog} GKlib m rt dl)
endforeach(prog)

set(CMAKE_C_COMPILER mpicc)
add_definitions(-DBDMPRUN)

if(BDMP_INSTALL)
  install(TARGETS bdmprun 
    RUNTIME DESTINATION bin)
endif()
