add_executable ( graphdiameter
      hll_counter.c
      parser.c
      graph.c
      hyper_anf_mpi.c
      parallel.c
)
target_link_libraries(graphdiameter bdmpi GKlib rt dl)
add_dependencies(graphdiameter ${GKlib_includes})

if(BDMP_INSTALL)
  install(TARGETS graphdiameter RUNTIME DESTINATION bin)
endif()
