add_executable ( serial
      hll_counter.c
      parser.c
      graph.c
      hyper_anf_serial.c
      serial.c
)

add_executable ( parallel
      hll_counter.c
      parser.c
      graph.c
      hyper_anf_mpi.c
      parallel.c
)

add_executable ( mpiHello
      mpi_hello.c
)

foreach(prog mpiHello
             parallel
             )
  target_link_libraries(${prog} bdmpi GKlib rt dl)
  add_dependencies(${prog} ${GKlib_includes})
endforeach(prog)
