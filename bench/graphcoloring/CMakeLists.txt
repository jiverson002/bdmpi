add_executable ( graphcoloring
      graphcoloring.c
)
target_link_libraries(graphcoloring bdmpi GKlib rt dl)
add_dependencies(graphcoloring ${GKlib_includes})

if(BDMP_INSTALL)
  install(TARGETS graphcoloring RUNTIME DESTINATION bin)
endif()
